import numpy as np
import pandas as pd
from ase.io import vasp

def gen_events(select_path, mode, inequival_sites=0):
    """
    生成存储所有事件的列表, 即一个位点下次跳跃发生时的可跳位置. 找到超胞中所有位点
    的所有可跳位点. 注意, 这些位点都是对称操作下等价的! 使用周期性边界条件.
    """
    # There have two modes for generate events list.
    # First way is generating events list by calculating the distance of two sites.
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
    # Third way is directly read from the csv file using pandas.
    elif mode == 'read':
        repeating_num_data = pd.read_csv("./repeating_num.csv")
        repeating_num = repeating_num_data.values.flatten().tolist()
        repeating_num = [int(value) for value in repeating_num]
        events_data = pd.read_csv("./events.csv")
        events = events_data.values.tolist()
        events = [[int(value) for value in row] for row in events]
    # 转换为ndarray格式。
    return repeating_num, np.array(events)