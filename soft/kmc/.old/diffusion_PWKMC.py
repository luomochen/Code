#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#-----------------------------------------
# 动力学蒙特卡洛模拟(KMC).
# 在无穷大晶格中所有等价位点下稀释极限中的扩散。
# 终止条件总共有两个, 分别是时间和步数。可根据
# 要求进行调整.
# 为了解决部分扩散路径的势垒较小，难以模拟长时间
# 扩散行为，通过Probability-Weighted Kinetic 
# Monte Carlo Method来解决这个问题。
#-----------------------------------------
import copy
import time
import numpy as np
import pandas as pd
import multiprocessing as mp
from numba import njit
from scipy.constants import physical_constants
from ase.io import vasp 

import parse

def list_extend(list, repeating_num):
    """
    存在部分等价的扩散路径, 需要对应的扩散势垒, 跳动几率等列表扩充.
    """
    list_ex = []
    for i in range(len(repeating_num)):
        list_ex.extend(repeating_num[i]*[list[i]])
    return list_ex

def gen_events(select_path, mode):
    """
    生成存储所有事件的列表, 即一个位点下次跳跃发生时的可跳位置. 找到超胞中所有位点
    的所有可跳位点. 注意, 这些位点都是对称操作下等价的! 使用周期性边界条件.
    """
    if mode == True:
        # 需要读取含有所有等价位点的超胞。
        sites = vasp.read_vasp("./sites.vasp")
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
                    if abs(path - dist) < 5E-4:
                        event.append(j)
                        num = num + 1
                repeating_num.append(num)
            path_num = sum(repeating_num)                
            event.extend([0]*path_num)
            events.append(event)
    # 将得到的所有event输出为csv文件, 便于检查.
        repeating_num_data = pd.DataFrame(repeating_num)
        repeating_num_data.to_csv("./repeating_num.csv", index=False)
        events_data = pd.DataFrame(events)
        events_data.to_csv("./events.csv", index=False)
    elif mode == 'read':
        repeating_num_data = pd.read_csv("./repeating_num.csv")
        repeating_num = repeating_num_data.values.flatten().tolist()
        repeating_num = [int(value) for value in repeating_num]
        events_data = pd.read_csv("./events.csv")
        events = events_data.values.tolist()
        events = [[int(value) for value in row] for row in events]
    # 转换为ndarray格式。
    return repeating_num, np.array(events)

def rate(T, barriers, jumpfreqs, pw, control_step):
    """
    根据DFT结果和Arrhenius公式计算指定温度下的反应速率和速率之和.
    使用float128保存数据.
    """
    total_rate = 0
    rates = []
    for i in range(len(barriers)):
        kb = physical_constants['Boltzmann constant in eV/K'][0]
        beta = 1.0/(kb * T)
        rate = np.float64(jumpfreqs[i] * np.exp(-barriers[i]*(beta)))
        rates.append(rate)
        total_rate = np.float64(total_rate + rate)
    # pw控制是否开启可能性加权动力学蒙特卡洛方法.
    if pw == True:
        weighted_factors = weighted_sampling(rates, total_rate, control_step)
        rates = rates / weighted_factors
        total_rate = np.sum(rates)
    elif pw == False:
        weighted_factors = np.ones(len(rates))
    name = ['rate', 'possibility', 'weighted_factors']
    rate_data = pd.DataFrame(columns=name, data=np.array([rates, np.float128(rates/total_rate), weighted_factors]).T)
    rate_data.to_csv(str(T)+'_rate.csv', index=False)
    return total_rate, np.array(rates)

def weighted_sampling(rates, total_rate, control_step):
    probilities = rates / total_rate
    weighted_factors = np.ones(len(probilities))
    # 将列表中重复值用set方法去除，并按照从大到小排序。获得第control_step个最大值。
    sorted_probilities = sorted(set(probilities), reverse=True)
    bigger_probilities = sorted_probilities[control_step]
    # 计算小于第control_step个最大值的概率之和。
    sum_prob = 0
    for prob in probilities:
        if prob < bigger_probilities:
            sum_prob += prob
    # 确定大于等于第control_step个最大值的概率的索引。
    indices = [index for index, prob in enumerate(probilities) if prob >= bigger_probilities]
    for index in indices:
        factor = round(probilities[index] / sum_prob)
        weighted_factors[index] = weighted_factors[index]*factor
    return weighted_factors

def cal_displacement(events, events_num):
    """
    后处理步骤, 计算累计位移, 方均根位移和扩散系数.
    """
    sites = vasp.read_vasp("./sites.vasp")
    displacement = np.zeros(3)
    for event in events:
        for i in range(events_num):
            # 每个事件下位移矢量和发生次数相乘获得该事件的累积位移.
            distance = sites.get_distance(event[0],event[i+1],mic=True,vector=True)
            distance = distance*event[i+1+events_num]
            displacement = displacement + distance
    return displacement

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
    return t, i, events

def kmc_loop(result_queue, t_control, nsteps, events, events_num, ini, total_rate, rates) -> None:
    """
    为了避免numba和多线程冲突, 再定义一个函数将kmc主循环包装起来. 先编译kmc_loop再
    将结果传递给队列.
    """
    t, i, events= kmc_iteration(t_control, nsteps, events, events_num, ini, total_rate, rates)
    result_queue.put([t, i, events])

def main_loop(pw, control_step, t_control, nsteps, repeat_run, T, repeating_number, events_origin, events_num, barriers, jumpfreqs) -> None:
    print(f"Start {T} K main loop.")
    # 计算运行时间.
    start=time.time()
    # 初始化输入.
    step_outcome = []
    barriers_ex = list_extend(barriers, repeating_number)
    jumpfreqs_ex = list_extend(jumpfreqs, repeating_number)
    total_rate, rates = rate(T, barriers_ex, jumpfreqs_ex, pw, control_step)
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
        #events_data = pd.DataFrame(result[2])
        #events_data.to_csv("./events_new.csv", index=False)
        displacement = cal_displacement(result[2], events_num)
        step_outcome.append([result[0], result[1], displacement[0], displacement[1], displacement[2]])
    end = time.time()
    print(f"{T} K finished, spend time: {round(end-start, 5)} s.")
    name = ['t', 'steps', 'd_x', 'd_y', 'd_z']
    outcome = pd.DataFrame(columns=name, data=step_outcome)
    outcome.to_csv(str(T)+'.csv', index=False)

def main():
    # 输入数据.
    pw, events_mode, control_step, repeat_run, nsteps, quantum_correction, \
            T_list, barriers, select_path = parse.parse_input('input.yaml')
    t_control_list = np.ones(len(T_list))*1E30
    # 初始化事件.
    repeating_number, events_origin = gen_events(select_path, event_gen_mode)
    events_num = sum(repeating_number)
    # 预先调用kmc_iteration保证numba提前对kmc和kmc_iteration编译, 确保不和多线程产生冲突.
    _ = kmc_iteration(1, 1, np.array([[0, 1, 0],[1, 0, 0]]), 1, 0, 1E12, [1E12])
    for i in range(len(T_list)):
        main_loop(pw, control_step, t_control_list[i], nsteps, repeat_run, T_list[i], repeating_number, events_origin, events_num, barriers, jumpfreqs)

if __name__ == "__main__":
    main()