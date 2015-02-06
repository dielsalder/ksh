import pdb as p
import split as s

mterf = p.Pdb('./pdb/3MVA.pdb')
dna = p.Dna(mterf.get_nucleic())
split = s.Split(dna)
split.iterbp(start = 12)
print split.fmin_bp()
