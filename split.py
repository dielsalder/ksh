import pdb as p
import dna
import ksh

pdb = p.Pdb('ideal.p.pdb')
pdb.get_all()
dna = p.Molecule(pdb)
