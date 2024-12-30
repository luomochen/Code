#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#-----------------------------------------
# 用于读取动力学蒙特卡洛模拟所需要的输入文件.
#-----------------------------------------
import yaml
import numpy as np

def read_input(filename='input.yaml'):
    with open(filename, "r") as file:
        input_data = yaml.safe_load(file)
    repeat_run = int(input_data["simulation_parameters"]["repeat_run"])
    nsteps = int(float(input_data["simulation_parameters"]["nsteps"]))
    T_list = input_data["simulation_parameters"]["temperatures"]
    barriers = input_data["simulation_parameters"]["diffusion_barriers"]
    jumpfreqs = np.float128(input_data["simulation_parameters"]["jump_frequencies"])*1E12
    select_path = input_data["simulation_parameters"]["distances"]
    return repeat_run, nsteps, T_list, barriers, jumpfreqs, select_path

def main() -> None:
    repeat_run, nsteps, T_list, barriers, jumpfreqs, select_path = read_input()
    print(repeat_run, nsteps, T_list, barriers, jumpfreqs, select_path)

if __name__ == "__main__":
    main()