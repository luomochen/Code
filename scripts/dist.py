#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#-----------------------------------------------
# Automatically setting KPOINTS and POTCAR
#-----------------------------------------------
import os
import re
from ase.io import read
import matplotlib.pyplot as plt
from ase.geometry import geometry

def get_pressure_from_dirname(dirname):
    """从目录名中提取压力值（如 '10GPa' -> 10.0）"""
    match = re.match(r"(\d+\.?\d*)GPa", dirname)
    if match:
        return float(match.group(1))
    return None

def get_shortest_oo_distance(atoms, distance_range=None):
    """计算结构中 O–O 原子对的最短距离，并按范围过滤"""
    o_indices = [i for i, atom in enumerate(atoms) if atom.symbol == "O"]
    min_dist = float("inf")
    for i in range(len(o_indices)):
        for j in range(i + 1, len(o_indices)):
            dist = atoms.get_distance(o_indices[i], o_indices[j], mic=True)
            if dist < min_dist:
                min_dist = dist
    if distance_range is not None:
        if not (distance_range[0] <= min_dist <= distance_range[1]):
            return None
    return min_dist

def collect_data_from_directories(distance_range=None):
    """从目录收集O–O最短距离，添加距离范围限制"""
    data = []
    for dirname in os.listdir("."):
        if os.path.isdir(dirname) and "GPa" in dirname:
            pressure = get_pressure_from_dirname(dirname)
            contcar_path = os.path.join(dirname, "CONTCAR")
            if pressure is not None and os.path.exists(contcar_path):
                try:
                    atoms = read(contcar_path, format="vasp")
                    oo_dist = get_shortest_oo_distance(atoms, distance_range)
                    if oo_dist is not None:
                        data.append((pressure, oo_dist))
                        print(f"{dirname}: {pressure} GPa -> O–O = {oo_dist:.3f} Å")
                    else:
                        print(f"{dirname}: O–O 距离不在范围内，跳过")
                except Exception as e:
                    print(f"读取 {dirname} 时出错: {e}")
    return sorted(data)

def plot_oo_distance_vs_pressure(data):
    """绘制 O–O 最短距离 vs 压力 的曲线图"""
    pressures, distances = zip(*data)

    plt.figure(figsize=(8, 6))
    plt.plot(pressures, distances, marker='o', linestyle='-')
    plt.xlabel("Pressure (GPa)")
    plt.ylabel("Shortest O–O distance (Å)")
    plt.title("O–O Shortest Distance vs. Pressure")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("oo_distance_vs_pressure.png", dpi=300)
    plt.show()

def main():
    distance_range = (3, 3.5)  # 用户可以在此处修改距离范围
    data = collect_data_from_directories(distance_range)
    if data:
        plot_oo_distance_vs_pressure(data)
    else:
        print("没有找到符合距离范围的有效数据。")

if __name__ == "__main__":
    main()