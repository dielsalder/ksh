import pdb
import ksh

class Fitter:
    def __init__(self, filename):
        self.pdb = pdb.Pdb(filename)
        self.atoms = self.pdb.get_nucleic()
        self.molecule = pdb.Dna(pdb.Molecule(self.atoms))

    def fit_selection(self, s):
        self.selection = self.molecule.select(s)
        self.selection.minimize()

