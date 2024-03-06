from ase.io import write,read
import numpy as np
import pymatgen.core as mg
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
import os
import argparse
from itertools import combinations
"""
description: 检测多个结构是否为等价，支持多种文件格式
usage:
python *.py POSCAR1 POSCAR2 POSCAR3 # POSCAR*为结构文件的路径
"""
#用法：python StruEqui.py POSCAR_*
def _ase2pymatgen(ase_stru):
    """将ase的结构转换为pymatgen的结构"""
    write(filename=".vasp", images=ase_stru, format="vasp")
    pymatgen_stru = Structure.from_file(filename=".vasp")
    os.remove(".vasp")
    return pymatgen_stru

def symmetryCheck(atoms1, atoms2):
    
    atoms1 = _ase2pymatgen(read(atoms1))
    atoms2 = _ase2pymatgen(read(atoms2))
   
    # primitive_cell=True： 可以比较不同原子数的两个结构；scale：是否缩放到同一体积，高精度下关闭此项
    comp = StructureMatcher(ltol=0.1, stol=0.1, angle_tol=1, primitive_cell=True, scale=False)
    return comp.fit(atoms1, atoms2)

parser = argparse.ArgumentParser(description="None")
parser.add_argument('POSCAR', nargs="*", help="Specify the POSCAR file")
prm = parser.parse_args()

dict_ = {}  # keys: 被检测的结构组合， values: 是否为相同结构,bool类型
for cb in list(combinations(prm.POSCAR,2)):
    results = symmetryCheck(cb[0], cb[1])
    dict_[cb] = results
output_file = "output.txt"
if np.array(list(dict_.values())).all():
    POSCAR_names = str(prm.POSCAR).strip("[]")
    with open(output_file, "w") as f:
        f.write("%s 均为相同结构\n" % POSCAR_names)
else:
    with open(output_file, "w") as f:
        f.write(">>> 存在非相同结构 <<<\n")
    for key, value in dict_.items():
        with open(output_file, "a") as f:
            f.write("%s 和 %s 是否为相同结构：%s\n" % (key[0], key[1], value))
#grep True output.txt >> compare
#gawk '{print $3}' compare >>data 
#if np.array(list(dict_.values())).all():
#    POSCAR_names = str(prm.POSCAR).strip("[]")
#    print("%s 均为相同结构" % POSCAR_names)
#else:
#    print(">>> 存在非相同结构 <<<")
#    for key,value in dict_.items():
#    
#    print("%s 和 %s 是否为相同结构：%s" % (key[0], key[1], value))
# ... 其他部分不变 ...
# ... 其他部分不变 ...
# ... 其他部分不变 ...

# ... 其他部分不变 ...

# 输出文件的路径
#output_file = "output.txt"
#
#if np.array(list(dict_.values())).all():
#    POSCAR_names = str(prm.POSCAR).strip("[]")
#    with open(output_file, "w") as f:
#        f.write("%s 均为相同结构\n" % POSCAR_names)
#else:
#    with open(output_file, "w") as f:
#        f.write(">>> 存在非相同结构 <<<\n")
#    for key, value in dict_.items():
#        with open(output_file, "a") as f:
#            f.write("%s 和 %s 是否为相同结构：%s\n" % (key[0], key[1], value)




