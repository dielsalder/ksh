# new version
import ksh
import base_pairs_i as bp
import base_pairs as bpn
import numpy as np
import sys

def get_crds(crd_file):
"""Solve from file of whitespace delimited xyz coordinates"""
return ksh.crdset(np.loadtxt(crd_file))

def parse_pdb(pdb_file):
"""Solve fit from pdb file"""
crdset = ksh.crdset(bp.strip_crds_strict_pdb(pdb_file))
return crdset

class TestCrds:
    """Coordinates to test"""
    def __init__(self, inpcrds, type='crds'):
        if type == 'crds':
	    self.crdset = get_crds(name)
	elif type == 'pdb':
	    self.crdset = parse_pdb(inpcrds)
	else:
            print "File type is not appropriate"

    def test_crds(self, title="Results"):
        """Fit existing crdset without minimizing"""
        self.crdset.calc_all()
        self.crdset.print_results(title)

    def fit_rotation(self, title="Best fit"):
        """Test minimization on crdset"""
        best = ksh.best_rotation(self.crdset)
        best.print_results(title)

    def test_rotation(self, reference_crds=self.crdset):
        """Print known fit then minimization attempt"""
        reference_crds.test_crds(title="True best coordinates")
        self.crdset.rotation_fit(title="Best fit")

#a = parse_pdb('bdna.4bp.bbatoms.pdb')
#test_crds(a, title="4 base pairs, backbone atoms")
#
#c = get_crds('bdna.pdb')
#test_crds(c)
#
#btest = get_crds('B.test.crd')
#test_crds(btest, title="One strand")
#rotation_fit(btest)

e55 = TestCrds('e55dna.pdb', type='pdb')
e55.test_crds

#idealb = parse_pdb('ideal_bdna.pdb')
#rotation_fit(idealb, title="Best fit from pdb")

#bp_idealb = bp.dna('ideal_bdna.pdb')
#bp_idealb.strip_strand_crds('id_bdna_strand.crd')
#idb_strand = get_crds('id_bdna_strand.crd')
#test_crds(idb_strand)
