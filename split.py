from numpy import fmin_bfgs

import pdb
import dna
import ksh

class Split(pdb.Dna):
    def __init__(self):
        self.all_res = {}

    def split(self, ibp):
        """Split helix at base pair ibp"""
        self.split = ibp
        self.front = pdb.crdset([a for a in self.atoms if a['pair'] < ibp])
        self.rear = pdb.crdset([a for a in self.atoms if a['pair'] > ibp])
        return self.front, self.rear

    def eval_res(self):
        """Calculate fit for current split and store res in self.all_res"""
        res_f = self.minimize((self.front).res)
        res_r = self.minimize((self.rear).res)
        res_sum = res_f + res_r
        self.all_res[self.split] = sum_res
        return sum_res

    def iterbp(self):
        """Fill all_res with res from each bp and find min"""
        for ibp in (2, self.numbp - 2):
            self.split(ibp)
            self.eval_res()
        self.res_min = min(self.all_res, key = self.all_res.get)
