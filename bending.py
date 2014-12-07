import base_pairs
import ksh

def helix_segments(base_pairs):
    # length of individual axis
    n = 5
    for pair in base_pairs:
        while base_pairs.index(pair) < (len(base_pairs)-n):
            helix = base_pairs[base_pairs.index(pair):base_pairs.index(pair)+n]
            crds =(helix[0] + helix[1])[2:]
            print crds
#           crds = (helix[0] + helix[1])[3:][0]
            yield ksh.crdset(crds)

def res_search(inp_pairs):
    lowest_1 = helix_segments(inp_pairs)
    lowest_2 = helix_segments(inp_pairs)
    for helix in helix_segments(inp_pairs):
        if helix.res < lowest_1.res:
            lowest_1 = helix
        elif helix.res < lowest_2.res:
            lowest_2 = helix
    return [lowest_1, lowest_2]

def add_pairs():
    pass

def find_best_pitch(inp_dna):
    print res_search(inp_dna)[0].res
    for min_res in res_search(inp_dna):
        add_pairs(min_res)

inp_pdb = './pdb/3MVA.pdb'
inp_dna = base_pairs.dna(inp_pdb)
find_best_pitch(inp_dna.base_pairs)
