#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#-----------------------------------------
# 动力学蒙特卡洛模拟(kMC).
# 在无穷大晶格中所有等价位点下稀释极限中的扩散。
# 终止条件总共有两个, 分别是时间和步数。可根据
# 要求进行调整.
#-----------------------------------------
import copy
import math
import time
import numpy as np
import pandas as pd
import multiprocessing as mp
from numba import njit
from scipy.constants import physical_constants
from ase.io import vasp 

def list_extend(list, repeating_num):
    """
    存在部分等价的扩散路径, 需要对应的扩散势垒, 跳动几率等列表扩充.
    """
    list_ex = []
    for i in range(len(repeating_num)):
        list_ex.extend(repeating_num[i]*[list[i]])
    return list_ex

def gen_events(select_path, output):
    """
    生成存储所有事件的列表, 即一个位点下次跳跃发生时的可跳位置. 找到超胞中所有位点
    的所有可跳位点. 注意, 这些位点都是对称操作下等价的! 使用周期性边界条件.
    """
    # 需要读取含有所有等价位点的超胞。
    sites = vasp.read_vasp("./Hi-_sup.vasp")
    site_number = sum(sites.numbers)
    # events的存储格式。[[],[],[]]，第一列是当前位点，后面的n列可跳位点。n列之后是
    # 对应的跳跃次数。为了保证能够使用numba加速，使用ndarray格式。
    events = []
    for i in range(site_number):
        event = [i]
        repeating_num = []
        for path in select_path:
            num = 0
            for j in range(site_number):
                dist = sites.get_distance(i, j, mic=True)
                # 判定等价路径的条件, 距离相等.         
                if abs(path - dist) < 1E-3:
                    event.append(j)
                    num = num + 1
            repeating_num.append(num)
        path_num = sum(repeating_num)                
        event.extend([0]*path_num)
        events.append(event)
    # 将得到的所有event输出为csv文件, 便于检查.
    if output == True:
        events_data = pd.DataFrame(events)
        events_data.to_csv("./events.csv", index=False)
    # 转换为ndarray格式。
    return repeating_num, np.array(events)

def rate(T, barriers, jumpfreqs):
    """
    根据DFT结果和Arrhenius公式计算指定温度下的反应速率和速率之和.
    使用float64保存数据.
    """
    total_rate = 0
    rates = []
    for i in range(len(barriers)):
        kb = physical_constants['Boltzmann constant in eV/K'][0]
        beta = 1.0/(kb * T)
        rate = np.float64(jumpfreqs[i] * np.exp(-barriers[i]*(beta)))
        rates.append(rate)
        total_rate = np.float64(total_rate + rate)
    return total_rate, np.array(rates)

def cal_diffusivity(events, events_num, t):
    """
    后处理步骤, 计算累计位移, 方均根位移和扩散系数.
    """
    sites = vasp.read_vasp("./Hi-_sup.vasp")
    displacement = np.zeros(3)
    for event in events:
        for i in range(events_num):
            # 每个事件下位移矢量和发生次数相乘获得该事件的累积位移.
            distance = sites.get_distance(event[0],event[i+1],mic=True,vector=True)
            distance = distance*event[i+1+events_num]
            displacement = displacement + distance
    # 三维扩散, 计算三个方向上的扩散系数再取平均.
    x_square = np.float64(math.pow(displacement[0], 2) * math.pow(10, -20))
    y_square = np.float64(math.pow(displacement[1], 2) * math.pow(10, -20))
    z_square = np.float64(math.pow(displacement[2], 2) * math.pow(10, -20))
    diff_xyz = np.array([x_square, y_square, z_square])/(2*t)
    diff = np.sum(diff_xyz)/3
    return diff, diff_xyz, displacement

@njit
# 和kmc相关的代码使用机器码.
def kmc(events, events_num, ini, total_rate, rates):
    """
    kmc过程主代码。
    """
    # 1. 随机安放一个H原子(初始化步骤), 该步骤位于主函数中.
    event = events[ini]
    # 2. 在0-1中获取一个随机数. 和累积函数相乘（R=sum(r_i)）,确定
    # rho1*R所在的区间, 确定所要发生的跳跃事件.
    rho1 = np.random.random()
    target = rho1 * total_rate
    sum_l = rates[0]
    if target < sum_l:
        fin = event[1]
        # 注意, 此时的rates列表已经完成扩充, 长度对应的所有事件数目.
        event[1+events_num] = event[1+events_num] + 1
    else:
        for i in range(2, events_num+1):
            sum_l = sum_l + rates[i-1]
            if target < sum_l:
                fin = event[i]
                event[i+events_num] = event[i+events_num] + 1
                break
    events[ini] = event
    # 3. 在0-1中获取一个随机数, 计算反应时间.
    rho2 = np.random.random()
    dt = -np.log(rho2) / total_rate
    #print(rho2, dt)
    # 返回更新后的时间, 位置和事件计数.
    return dt, fin, events

@njit
def kmc_iteration(t_control, nsteps, events, events_num, ini, total_rate, rates):
    """
    kmc过程主循环. 循环终止条件可选择抵达累积的物理时间和抵达累积的步数.
    """
    # t_control控制累积时间, nsteps确定最少循环次数.
    t = 0
    i = 0
    while i < nsteps:
        dt, ini, events = kmc(events, events_num, ini, total_rate, rates)
        t = t + dt
        i = i + 1
        if t > t_control:
            break
    return t, events, i

def kmc_loop(result_queue, t_control, nsteps, events, events_num, ini, total_rate, rates):
    """
    为了避免numba和多线程冲突, 再定义一个函数将kmc主循环包装起来. 先编译kmc_loop再
    将结果传递给队列.
    """
    t, events, i= kmc_iteration(t_control, nsteps, events, events_num, ini, total_rate, rates)
    result_queue.put([t, events, i])

def main_loop(t_control, nsteps, repeat_run, T, repeating_number, events_origin, events_num, barriers, jumpfreqs):
    # 计算运行时间.
    start=time.time()
    # 初始化输入.
    t_list = []
    diff_list = []
    step_outcome = []
    barriers_ex = list_extend(barriers, repeating_number)
    jumpfreqs_ex = list_extend(jumpfreqs, repeating_number)
    total_rate, rates = rate(T, barriers_ex, jumpfreqs_ex)
    # 在kMC循环处引入多线程, 即repeat_run次重复分别放在不同的处理器上跑.
    # 默认进程数为调用的核数.
    manager = mp.Manager()
    result_queue = manager.Queue()
    pool = mp.Pool()
    for run in range(repeat_run):
        events = copy.deepcopy(events_origin)
        # 随机在一个位点上生成一个H原子.
        ini = np.random.randint(0, len(events_origin))
        pool.apply_async(func=kmc_loop, args=(result_queue, t_control, nsteps, events, events_num, ini, total_rate, rates))
    pool.close()
    pool.join()
    # 放入终止标志.
    result_queue.put(None)
    results = []
    # 获取结果并处理结果.
    while True:
        result = result_queue.get()
        if result is None:
            break
        results.append(result)
    for result in results:
        diff, diff_xyz, displacement = cal_diffusivity(result[1], events_num, result[0])
        t_list.append(result[0])
        diff_list.append(diff)
        step_outcome.append([result[0], diff,
                             diff_xyz[0], diff_xyz[1], diff_xyz[2], 
                             displacement[0], displacement[1], displacement[2],
                             result[2]])

    diff_mean = np.average(np.array(diff_list), weights=np.array(t_list))
    diff_std = np.std(np.array(diff_list))
    end = time.time()
    print(f"{T}: mean diffusivity: {diff_mean}, SD: {diff_std}, spend time: {end-start} s.")
    name = ['t', 'D', 'D_x', 'D_y', 'D_z', 'd_x', 'd_y', 'd_z', 'steps']
    outcome = pd.DataFrame(columns=name, data=step_outcome)
    outcome.to_csv(str(T)+'.csv', index=False)
    return [T, diff_mean, diff_std]

def main():
    # 输入数据.
    T_list = [300, 500, 700, 1000, 1500, 2000, 2500, 3000]
    repeat_run = int(20000)
    #t_control_list = np.ones(8)*1E8
    t_control_list = np.float64([3.654, 8.89E-3, 6.74E-4, 9.74E-5, 2.165E-5, 1.019E-5, 6.492E-6, 4.781E-6])
    nsteps = int(1E7)
    select_path = [3.10144, 3.10266, 3.16324]
    barriers = [2.61355, 0.38895, 2.62301]
    jumpfreqs = np.float64([74.932E12, 4.67931E12, 74.22347E12])
    # 初始化事件.
    repeating_number, events_origin = gen_events(select_path, True)
    events_num = sum(repeating_number)
    # 预先调用kmc_iteration保证numba提前对kmc和kmc_iteration编译, 确保不和多线程产生冲突.
    _ = kmc_iteration(1, 1, np.array([[0, 1, 0],[1, 0, 0]]), 1, 0, 1E12, [1E12])
    diff_list = []
    for i in range(len(T_list)):
        outcome = main_loop(t_control_list[i], nsteps, repeat_run, T_list[i], repeating_number, events_origin, events_num, barriers, jumpfreqs)
        diff_list.append(outcome)
    name = ['T', 'mean diffusivity', 'SD']
    outcome = pd.DataFrame(columns=name, data=diff_list)
    outcome.to_csv('outcome.csv', index=False)

if __name__ == "__main__":
    main()