import sys
import numpy
import base_pairs_i as bp
import ksh

def strip_crds_strict(atom):
    str_crds = atom[2:5]
    return [float(x) for x in str_crds]

def get_pair_crds(dna_mol, num):
    crds = []
    for pair in dna_mol.base_pairs[0:num+1]:
        for nuc in pair:
            for atom in nuc:
                crds.append(strip_crds_strict(atom))
    return crds

inp_pdb = './pdb/ideal_bdna.pdb'
dna_mol = bp.dna(inp_pdb)

n_res = []
for n in range(1, 9):
    crds = numpy.array(get_pair_crds(dna_mol, n))
    n_crdset = ksh.crdset(crds)
    n_res.append([n, n_crdset.res])

#all_pairs_res = ksh.crdset(dna_mol.all_crds).res

#for n in n_res:
#    print n[1]
#dna_crdset = ksh.crdset(dna_mol.all_crds)
#dna_crdset.calc_all
#ksh.print_results(dna_crdset, "DNA molecule")
#print numpy.array(n_res)
