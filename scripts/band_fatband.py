#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#---------------------------------------------
# A script for ploting the fat band structure.
#---------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from pyprocar import bandsplot
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.plotter import BSPlotter, BSPlotterProjected

class band():
    """docstring for band."""
    def __init__(self, filename) -> None:
        self.filename = filename

    def get_band_data(self):
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

    def n_band_eig_range(self, bs_data):
        branch_number = len(bs_data['distances'])
        band_number = len(bs_data['energy']['1'][0])
        eig_ranges = []
        for i in range(band_number):
            branch_extremums = []
            for j in range(branch_number):
                branch_max = max(bs_data['energy']['1'][j][i])
                branch_min = min(bs_data['energy']['1'][j][i])
                branch_extremums.append(branch_max)
                branch_extremums.append(branch_min)
            eig_ranges.append([min(branch_extremums), max(branch_extremums)])
        with open ('band_eig_value_range.dat', 'w') as file:
            for i in range(band_number):
                file.write(f'band{i+1}: {eig_ranges[i]}\n')
        return eig_ranges
    
    def band_structure_plot(self, bs_data, is_spin_polarized):
        # 设置颜色
        color_up = "blue"  # 自旋向上的颜色.
        color_down = "red"  # 自旋向下的颜色.
        # 绘制能带图.
        fig, ax = plt.subplots(figsize=(10, 8), dpi=150)
        # 绘制每一条能带.
        branch_number = len(bs_data['distances'])
        band_number = len(bs_data['energy']['1'][0])
        #colors = plt.cm.jet(np.linspace(0, 1, band_number))
        # 如果是自旋极化，绘制自旋向下的带.
        if is_spin_polarized:
            for i in range(branch_number):  # '1' 表示自旋向上的带.
                for j in range(band_number):
                    ax.plot(bs_data['distances'][i], 
                            bs_data['energy']['1'][i][j], 
                            color=color_up, 
                            label='Spin Up')
            for i in range(branch_number):  
                for j in range(band_number):
                    ax.plot(bs_data['distances'][i], 
                            bs_data['energy']['-1'][i][j], # '-1' 表示自旋向下的带.
                            color=color_down, 
                            label='Spin Up')
        else:
            for i in range(branch_number):
                for j in range(band_number):
                    if j > 7:
                        color = color_up
                    else:
                        color = color_down
                    if j < 11:
                        ax.plot(bs_data['distances'][i], 
                                bs_data['energy']['1'][i][j], 
                                color=color,
                                #label='Band'+str(j)
                                )
        # 设置x轴.
        ax.set_xticks(bs_data['ticks']['distance'])
        ax.set_xticklabels(bs_data['ticks']['label'])
        ax.set_xlim(bs_data['ticks']['distance'][0], bs_data['ticks']['distance'][-1])
        ax.set_ylim(-7, 22)
        # 设置y轴.
        ax.set_ylabel("Energy (eV)")
        # 添加网格.
        ax.grid(True)
        # 添加图例.
        ax.legend()
        fig.savefig('./band.png')
        plt.show()

def fatband(element_list, orbital_list):
    vasprun = Vasprun('vasprun.xml', parse_projected_eigen=True)
    band_structure = vasprun.get_band_structure(line_mode=True)
    projected_plotter = BSPlotterProjected(band_structure)
    for element in element_list:
        for orbital in orbital_list:
            ax = projected_plotter.get_projected_plots_dots({element:[orbital]}, 
                                                            zero_to_efermi=False,
                                                            ylim=(-10, 18))
            plt.savefig(element+orbital+'.png')
            #plt.show()

def pyprocar_fatband():
    orbitals = ['s', 'px', 'py', 'pz']
    elements = ['Si', 'Sb']
    atoms = [[0, 1], [2, 3, 4, 5]]
    for i in range(len(atoms)):
        for j in range(len(orbitals)):
            bandsplot(code='vasp', dirname='./', fermi=-2.44727085,
                                mode='parametric', elimit=[-5, 2], 
                                kticks=[0,59,119,179], knames=['G','K','M','G'],
                                cmap='jet', orbitals=[j], atoms=atoms[i],
                                savefig=elements[i]+orbitals[j]+'.png')

def main():
    vasprun_file = band('vasprun.xml')
    bs_data, is_spin_polarized = vasprun_file.get_band_data()
    vasprun_file.n_band_eig_range(bs_data)
    vasprun_file.band_structure_plot(bs_data, is_spin_polarized)
    #element_list = ['Mo', 'N', 'Si']
    #element_list = ['Si']
    #orbital_list = ['s', 'p', 'd']
    #fatband(element_list, orbital_list)

if __name__ == "__main__":
    main()
