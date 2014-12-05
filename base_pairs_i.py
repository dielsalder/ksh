import re
import numpy as np

def strip_crds(pdb):
    """Strip coordinates, atom names and residue numbers"""
    crds = []
    for line in pdb:
        crds.append(line[2:4] + [line[5]] + line[6:10])
    return crds

def strip_crds_strict_pdb(pdb_file):
    """Extract only xyz coordinates from a pdb file"""
    with open(pdb_file) as pdb:
        crds = [line.split()[5:8] for line in pdb if line.startswith('ATOM')]
        crds = [[float(x) for x in atom] for atom in crds]
    return crds

class dna:
    # extract coordinates from input pdb file
    def get_atoms(self, inp_pdb):
        self.atoms = []
        with open(inp_pdb) as pdb:
            for line in pdb:
                if re.search('^ATOM', line):
                    row = line.split()
                    self.atoms.append(row)
        return self.atoms

    def mask_atoms(self):
#       self.masked_atoms = [atom if atom[2] in self.mask for atom in self.atoms]
#       return self.masked_atoms
        pass

    def group_nucs(self, crds):
        prev_a = '        '
        res = []
        nucs = []
        for a in crds:
            if a[1] == prev_a[1]:
                res.append(a)
            else:
                nucs.append(res)
                res = []
            prev_a = a
        return nucs

    def parse_strand(self):
        """Store coordinates of strands separately"""
        prev_line = '11111111'
        for line in self.atoms:
            # [4] is the residue number
            if line[4] != prev_line[4] and line[4] == '1':
                i = self.atoms.index(line)
                self.strand1 = self.atoms[0:i]
                self.strand2 = self.atoms[i:]
            prev_line = line
        return [self.strand1, self.strand2]

    def pair_crds(self):
        """Create list of pairs of nucleotides"""
        self.base_pairs = []
        self.crds_a = strip_crds(self.strand1)
        self.crds_b = strip_crds(reversed(self.strand2))
        self.nucs_a = self.group_nucs(self.crds_a)
        self.nucs_b = self.group_nucs(self.crds_b)
        for (a, b) in zip(self.nucs_a, self.nucs_b):
            pair = [a, b]
            self.base_pairs.append(pair)
        return self.base_pairs

    def print_strands(self):
        for (a, b) in zip(self.strand1, reversed(self.strand2)):
            if a[0] == 'C1' or a[0] == 'C1\'':
                print a[5], a[3], '--', b[3], b[5]

    def strip_strand_crds(self, out):
        """Print coordinates for individual strands in separate files"""
        with open(out, 'w') as out:
            for atom in self.strand1:
                out.write(' '.join(atom[5:8]))
                out.write("\n")

    def printf_pairs(self, out):
        self.pair_crds()
        with open(out, 'w') as out:
            for pair in self.base_pairs:
                nuc_a = pair[0]
                nuc_b = pair[1]
                for (a, b) in zip(nuc_a, nuc_b):
                    out.write(('{} {} {} {} {} {}    {} {} {} {} {} {}\n').format(a[0], a[1], a[2], a[3], a[4], a[5], b[0], b[1], b[2], b[3], b[4], b[5]))
                out.write('\n')

    def __init__(self, inp_pdb):
        self.get_atoms(inp_pdb)
        self.mask = ('C3\'', 'O4\'', 'C5\'', 'O3\'', 'P', '05\'')
       #self.atoms = self.mask_atoms()
        self.parse_strand()
        self.pair_crds()

inp_pdb = 'ideal_bdna.pdb'
out = 'bdna_pairs.dat'
#a = dna(inp_pdb)
#a.printf_pairs(out)
