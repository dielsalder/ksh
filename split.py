import pdb as p
import dna
import ksh

class Split(pdb.Dna):
    def __init__(self):
        self.s_fits = {}

    def split(self, i):
        front = [a for a in self.atoms if a['pair'] < i]
        rear = [a for a in self.atoms if a['pair'] > i]
        yield front, rear

    def gen_split(self, np):
        for i in range(1, np):
            yield self.split(i)

    def fit_split(self, segs):
        for s in segs:
            p = s[0]['pair']
            s = s[0]['strand']
            fit = self.minimize(s)
            self.s_fits.append([p, s] : fit)

def test(filename):
    pdb = p.Pdb(filename)
    dna = p.Dna(p.Molecule(pdb.get_all))

