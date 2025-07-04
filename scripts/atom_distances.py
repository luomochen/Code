#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#---------------------------------------------
# A script for ploting the band structure.
#---------------------------------------------
from ase.io import read
import numpy as np
# 读取 POSCAR 文件
structure = read("CONTCAR", format="vasp")

# 提取所有氢原子和非氢原子的索引
h_indices = [i for i, atom in enumerate(structure) if atom.symbol == 'H']
non_h_indices = [i for i, atom in enumerate(structure) if atom.symbol != 'H']

for h_idx in h_indices:
    # 计算当前氢原子到所有非氢原子的距离
    distances = structure.get_distances(h_idx, non_h_indices, mic=True)
    
    # 按距离排序并取前五
    sorted_indices = np.argsort(distances)
    top5_indices = sorted_indices[:5]
    top5_distances = distances[top5_indices]
    
    print(f"氢原子 {h_idx} 的前五近邻（仅非氢原子）:")
    for rank, (list_idx, dist) in enumerate(zip(top5_indices, top5_distances), 1):
        original_idx = non_h_indices[list_idx]  # 映射回原结构的原子索引
        symbol = structure[original_idx].symbol
        print(f"  第{rank}近: {symbol} (索引 {original_idx}), 距离 {dist:.3f} Å")