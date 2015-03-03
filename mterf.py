import pdb as p
import split as s

mterf = p.Pdb('./pdb/3MVA.pdb')
dna = p.Dna(mterf.get_nucleic())
split = s.Split(dna)

#print "%3d  %12d\t%4d\t%4d\t%4d\t%4d" % (split.i_split, self.res_sum,
#            split.front.phi, self.front.the, self.rear.phi, self.rear.the)

split.iterbp()
print split.fmin_bp()
