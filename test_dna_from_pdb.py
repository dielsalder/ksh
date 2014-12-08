from dna_from_pdb import *

def print_pair_crds(dna):
    for crds in dna.pair_crds():
        for atom in crds:
            print atom
        print '\n'

#mterf = DnaMolecule('3MVA.pdb', mask=['C3\'', 'O4\'', 'C5\'', 'O3\'', 'P', '05\''])
bdna = DnaMolecule('e55dna.pdb', mask=['C3\'', 'O4\'', 'C5\'', 'O3\'', 'P', '05\''])
print_pair_crds(bdna)
