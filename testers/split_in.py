from pdb import *
from split import *
pdbfile = Pdb('./pdb/ideal_bdna.pdb')
mol = Molecule(pdbfile.get_all())
bdna = Dna(mol.atoms)

split = Split(bdna)
split.iterbp()
print split.fmin_bp()
print split.res_all
