import numpy as np
import ksh
import dna_from_pdb as dna

def segment(dna, num=1, max=22):
    pair_crds = []
    while num < max:
        pair_crds = pair_crds + (next(dna.pair_crds()))
        num += 1
	np_pair_crds = np.array(pair_crds)
        yield np_pair_crds, num

bdna = dna.DnaMolecule('idealbdna.bp.pdb')
#ksh.best_rotation(ksh.crdset(bdna.all_crds)).print_results()

print "\nnbp\tResidual\tRadius\t\tPitch"
for helix in segment(bdna):
    crds, num = helix
    crdset = ksh.best_rotation(ksh.crdset(crds))
#    best = ksh.best_rotation(crdset)
#    print num, "\t", best.res_v, "\t", best.radius, "\t\t", best.center
    crdset.calc_all()
    print num, "\t", crdset.res_v, "\t", crdset.radius, "\t", crdset.pitch
