#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import re
import yaml
import numpy as np

def get_calculated_point():
    """
    Identify folders matching 'stable', 'stable_*', 'saddle', or 'saddle_*'.
    These are considered calculation directories.
    """
    files = os.listdir()
    saddle_list = [f for f in files if re.match(r"^saddle(_\d+(\.\d+)?)?$", f)]
    stable_list = [f for f in files if re.match(r"^stable(_\d+(\.\d+)?)?$", f)]
    return stable_list, saddle_list

def get_frequency(folder, max_lines=100000):
    f_list, fi_list = [], []
    path = os.path.join(folder, "OUTCAR")
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            f.seek(0, os.SEEK_END)
            file_size = f.tell()
            block_size = 1024
            blocks, lines_found = [], 0
            while file_size > 0 and lines_found < max_lines:
                file_size = max(0, file_size - block_size)
                f.seek(file_size)
                block = f.read(block_size)
                blocks.append(block)
                lines_found = sum(block.count('\n') for block in blocks)
            tail = ''.join(reversed(blocks)).splitlines()[-max_lines:]
    except Exception:
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            tail = f.readlines()

    for line in tail:
        if "THz" in line and "f" in line:
            match_real = re.search(r"f\s+=\s*([-+]?[0-9]*\.?[0-9]+)", line)
            match_img = re.search(r"f/i=\s*([-+]?[0-9]*\.?[0-9]+)", line)
            if match_real:
                f_list.append(float(match_real.group(1)))
            if match_img:
                fi_list.append(float(match_img.group(1)))
    return np.array(f_list), np.array(fi_list)

def save_frequencies(name, f, fi):
    with open(f"{name}_f.yaml", "w", encoding="utf-8") as fout:
        yaml.dump(f.tolist(), fout, allow_unicode=True)
    with open(f"{name}_fi.yaml", "w", encoding="utf-8") as fio:
        yaml.dump(fi.tolist(), fio, allow_unicode=True)

def main():
    stable_dirs, saddle_dirs = get_calculated_point()
    for folder in stable_dirs + saddle_dirs:
        real_f, imag_f = get_frequency(folder)
        save_frequencies(folder, real_f, imag_f)
        print(f"Saved {folder}_f.yaml and {folder}_fi.yaml")

if __name__ == "__main__":
    main()