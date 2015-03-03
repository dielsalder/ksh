import scipy.optimize as sp
import pdb
import dna
import ksh

class Split(pdb.Dna):
    def __init__(self, Dna):
        self.atoms = Dna.atoms
        self.res_all = {}
        self.angles_all = {}
        self.set_numbp()

    def split(self, ibp):
        """Split helix at base pair ibp"""
        self.i_split = ibp
        self.front = pdb.Fit([a for a in self.atoms if a['pair'] < ibp])
        self.rear = pdb.Fit([a for a in self.atoms if a['pair'] > ibp])
        return self.front, self.rear

    def eval_res(self):
        """Calculate fit for current split and store res in self.res_all"""
        self.front.minimize()
        self.rear.minimize()
        res_f = self.front.res
        res_r = self.rear.res
        self.res_sum = res_f + res_r
        self.res_all[self.i_split] = self.res_sum

        self.angles_all[self.i_split] = ([[self.front.phi, self.front.the],
            [self.rear.phi, self.rear.the]])
        return self.res_sum

    def iterbp(self, start = 2):
        """Fill res_all with res from selected bp and print calculations"""
        print "i  res         \tfphi\tfthe\trphi\trthe" # header
        for ibp in range(start, self.numbp - 2):
            self.split(ibp)
            #self.eval_res()
            self.eval_res()
            print self
        return self.res_all

    def fmin_bp(self):
        """Find ibp of minimum residual in res_all"""
        pass
        self.bp_min = min(self.res_all, key = self.res_all.get)
        return self.bp_min

    def min_crdset(self):
        """Evaluate helix parameters for minimum split"""
        self.split(ibp)
        self.front.minimize()
        self.rear.minimize()
        return self.front, self.rear, self.i_split

    def __str__(self):
        return "%3d  %12d\t%4d\t%4d\t%4d\t%4d" % (self.i_split, self.res_sum,
                self.front.phi, self.front.the, self.rear.phi, self.rear.the)
