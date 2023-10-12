import sys
import ase.io.vasp

print(sys.argv[:])
x,y,z = [int(i) for i in sys.argv[1:4]]
cell = ase.io.vasp.read_vasp("POSCAR")
ase.io.vasp.write_vasp("POSCAR",cell*(x, y, z), direct=True,sort=True)
