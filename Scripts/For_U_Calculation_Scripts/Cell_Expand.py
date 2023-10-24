import sys
import ase.io.vasp

# 从stdin中读取的参数包含该脚本文件名，所以列表顺序应当为1:4而非0:4。
print(sys.argv[:])
x,y,z = [int(i) for i in sys.argv[1:4]]
cell = ase.io.vasp.read_vasp("POSCAR")
# 注意！此处的sort打开后POSCAR中的原子将按照元素周期表的顺序排列。
ase.io.vasp.write_vasp("POSCAR",cell*(x, y, z), direct=True,sort=True)
