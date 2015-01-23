from pdb import *
from split import *
pdbfile = Pdb('./pdb/ideal_bdna.pdb')
mol = Molecule(pdbfile.get_all())
bdna = Dna(mol.atoms)

split = Split(bdna)
split.split(5)
print split.eval_res()
# 2473.1115
