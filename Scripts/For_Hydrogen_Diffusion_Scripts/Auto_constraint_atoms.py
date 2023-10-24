import ase.io.vasp
from ase.constraints import FixAtoms

atoms = ase.io.vasp.read_vasp("CONTCAR")
c = FixAtoms(indices=[atom.index for atom in atoms])
atoms.set_constraint(c)
ase.io.vasp.write_vasp("POSCAR", atoms, direct=True)