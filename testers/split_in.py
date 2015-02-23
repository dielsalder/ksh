import pdb as p
import split as s

mterf = p.Pdb('./pdb/3MVA.pdb')
dna_all = p.Dna(mterf.get_nucleic())
dna_p = p.Dna(dna_all.select(atom_name = ['P']))

split = s.Split(dna_p)
split.iterbp()
print split.fmin_bp()
print split.res_all
