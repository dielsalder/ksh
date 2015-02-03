import scipy.optimize as sp
import pdb
import dna
import ksh

class Split(pdb.Dna):
    def __init__(self, Dna):
        self.atoms = Dna.atoms
        self.res_all = {}
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
        res_sum = res_f + res_r
        self.res_all[self.i_split] = res_sum
        return res_sum

    def iterbp(self):
        """Fill res_all with res from each bp"""
        for ibp in range(2, self.numbp - 2):
            self.split(ibp)
            self.eval_res()
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
        return self.res_sum
