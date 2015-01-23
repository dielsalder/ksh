from pdb import *
pdbfile = Pdb('./pdb/ideal_bdna.pdb')
mol = Molecule(pdbfile.get_all())
bdna = Dna(mol.atoms)
fit = Fit(bdna.atoms)
fit.minimize()
fit.crdset.calc_all()
fit.crdset.print_results()
