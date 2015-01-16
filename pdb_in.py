from pdb import *
pdbfile = Pdb('./pdb/ideal_bdna.pdb')
mol = Molecule(pdbfile.get_all())
bdna = Dna(mol.atoms)
#print [a for a in bdna.atoms if a['pair'] in range(1, 2)]
#print bdna.get_pairs(2, 3)

print bdna.get_pair(3)
