from numpy import fmin_bfgs

import pdb
import dna
import ksh

class Split(pdb.Dna):
    def __init__(self):
        self.res_all = {}

    def split(self, ibp):
        """Split helix at base pair ibp"""
        self.i_split = ibp
        self.front = pdb.crdset([a for a in self.atoms if a['pair'] < ibp])
        self.rear = pdb.crdset([a for a in self.atoms if a['pair'] > ibp])
        return self.front, self.rear

    def eval_res(self):
        """Calculate fit for current split and store res in self.res_all"""
        res_f = self.minimize((self.front).res)
        res_r = self.minimize((self.rear).res)
        res_sum = res_f + res_r
        self.res_all[self.i_split] = sum_res
        return sum_res

    def iterbp(self):
        """Fill res_all with res from each bp"""
        for ibp in (2, self.numbp - 2):
            self.split(ibp)
            self.eval_res()
        return self.res_all

    def fmin_res(self):
        """Find ibp of minimum residual in res_all"""
        self.bp_min = min(self.res_all, key = self.res_all.get)
        return self.bp_min

    def min_crdset(self):
        """Evaluate helix parameters for minimum split"""
        self.split(ibp)
        self.front.calc_all()
        self.rear.calc_all()
        return self.front, self.rear, self.i_split
