from pdb import *
pdbfile = Pdb('./pdb/ideal_bdna.pdb')
mol = Molecule(pdbfile.get_all())
bdna = Dna(mol.atoms)
