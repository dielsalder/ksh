import ksh
import base_pairs_i as bp
import numpy as np

def textcrds(crd_file):
    """Make crdset from file of whitespace delimited xyz coordinates"""
    return ksh.crdset(np.loadtxt(crd_file))

def pdbcrds(pdb_file):
    """Make crdset from pdb file"""
    crdset = ksh.crdset(bp.strip_crds_strict_pdb(pdb_file))
    return crdset

class Test:
    """Load and test set of coordinates"""
    def __init__(self, name, type='crds'):
        """Load coordinates from text or pdb"""
        if type == 'crds':
            self.crdset = textcrds(name)
        elif type == 'pdb':
            self.crdset = pdbcrds(name)
        else:
            print "Please specify appropriate file type"""

    def fit(self, title="Results"):
        """Fit existing crdset without minimizing"""
        self.crdset.calc_all()
        self.crdset.print_results(title)

    def min_fit(self, title="Best fit"):
        """Test minimization"""
        best = ksh.best_rotation(self.crdset)
        best.print_results(title)

   # def min_ref(self, reference_crds=self.crdset):
   #     """Print known fit then minimization attempt"""
   #     reference_crds.test_crds(title="True best coordinates")
   #     self.crdset.rotation_fit(title="Best fit")

   # def rotate(self, phi=180, the=0):
   #     """Fit rotation by specified phi and theta"""
   #     pass

bdna = Test('./pdb/bdna.pdb')
bdna.fit()
#while True:
#	cmd, *args = shlex.split(input('\> '))
