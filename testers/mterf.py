import pdb as p
import split as s

mterf = p.Pdb('./pdb/3MVA.pdb')
dna = p.Dna(mterf.get_nucleic())
split = s.Split(dna)
split.iterbp()
print split.fmin_bp()
