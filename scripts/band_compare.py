#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#---------------------------------------------
# A script for ploting the band structure.
#---------------------------------------------
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.plotter import BSPlotter

class band():
    """docstring for band."""
    def __init__(self, filename) -> None:
        self.filename = filename
        
    def get_dft_band_data(self):
        # 读取VASP输出.
        vasprun = Vasprun(self.filename)
        # 提取能带结构.
        band_structure = vasprun.get_band_structure(line_mode=True)
        # 检查是否包含自旋信息.
        is_spin_polarized = band_structure.is_spin_polarized
        # 创建能带图绘制对象.
        bs_plotter = BSPlotter(band_structure)
        bs_data = bs_plotter.bs_plot_data(zero_to_efermi=False)
        return bs_data, is_spin_polarized

    def dft_band_structure_plot(self, bs_data, is_spin_polarized):
        # 绘制每一条能带.
        branch_number = len(bs_data['distances'])
        band_number = len(bs_data['energy']['1'][0])
        fig, ax = plt.subplots(figsize=(10, 8), dpi=150)
        # 如果是自旋极化，绘制自旋向下的带.
        if is_spin_polarized:
            print("spin polarized!")
            for i in range(branch_number):  # '1' 表示自旋向上的带.
                for j in range(band_number):
                    ax.plot(bs_data['distances'][i], 
                            bs_data['energy']['1'][i][j], 
                            color='b', 
                            label='Spin Up')
            for i in range(branch_number):  # '-1' 表示自旋向下的带.
                for j in range(band_number):
                    ax.plot(bs_data['distances'][i], 
                            bs_data['energy']['-1'][i][j], 
                            color='r', 
                            label='Spin Down')
        else:
            for i in range(branch_number):
                for j in range(band_number):
                    ax.plot(bs_data['distances'][i], 
                            bs_data['energy']['1'][i][j],
                            color='b')
        # 设置x轴.
        ax.set_xticks(bs_data['ticks']['distance'])
        ax.set_xticklabels(bs_data['ticks']['label'])
        ax.set_xlim(bs_data['ticks']['distance'][0], bs_data['ticks']['distance'][-1])
        ax.set_ylim(-10, 15)
        # 设置y轴.
        ax.set_ylabel("Energy (eV)")
        # 添加网格.
        ax.grid(True)
        legend_lines = [Line2D([0], [0], color='b', linestyle='-', label='DFT band'),
                        Line2D([0], [0], color='r', linestyle='--', label='wannier band')]
        ax.legend(handles=legend_lines)
        return fig, ax
    
    def get_wannier_band_data(self):
        wannier_data = []  # 用于存储所有段的列表
        Kpoints = []
        energy = []
        with open(self.filename, 'r') as file:
            for line in file:
                stripped_line = line.strip()
                if not stripped_line:  # 空行表示段结束
                    if Kpoints:  # 如果当前段非空，存储并重置
                        wannier_data.append([Kpoints, energy])
                        Kpoints = []
                        energy = []
                else:
                    # 将当前行分割为浮点数并添加到当前段中
                    row = list(map(float, stripped_line.split()))
                    Kpoints.append(row[0])
                    energy.append(row[1])
        # 确保最后一段也被添加
        if Kpoints:
            wannier_data.append([Kpoints, energy])
        return wannier_data
    
    def wannier_band_structure_plot(self, fig, ax, wannier_data):
        for band in wannier_data:
            ax.plot(band[0], band[1], color='r', linestyle='--')
        return fig, ax 
               
    def save_fig():
        pass
               
def main():
    vasprun_file = band('vasprun.xml')
    wannier_file = band('wannier90_band.dat')
    bs_data, is_spin_polarized = vasprun_file.get_dft_band_data()
    fig, ax = vasprun_file.dft_band_structure_plot(bs_data, is_spin_polarized)
    wannier_data = wannier_file.get_wannier_band_data()
    fig, ax = wannier_file.wannier_band_structure_plot(fig, ax, wannier_data)
    fig.savefig('band.png', format='png', dpi=300, bbox_inches='tight')
    plt.show()
if __name__ == "__main__":
    main()