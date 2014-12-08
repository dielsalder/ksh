import re

def strip_crds(pdb):
    crds = []
    for line in pdb:
        crds.append(line[2:4] + [line[5]] + line[6:10])
    return crds

# only extract x, y and z
def strip_crds_strict(pdb):
    crds = []
    for line in pdb:
        crds.append(line[6:9])
    return crds

class dna:
    # extract coordinates from input pdb file
    def get_atoms(self, inp_pdb):
        self.atoms = []
        with open(inp_pdb) as pdb:
            for line in pdb:
                if re.search('(^ATOM)\s*\S*\s*\S*\s*(DA5|DA3|DA|DT5|DT3|DT|DG5|DG3|DG|DC5|DC3|DC)', line):
                    row = line.split()
                    self.atoms.append(row)
        return self.atoms

    def mask_atoms(self):
        self.masked_atoms = []
        for atom in self.atoms:
            if atom[2] in self.mask:
                self.masked_atoms.append(atom)
        return self.masked_atoms

    def group_nucs(self, crds):
        prev_a = ' 2      '
        res = []
        nucs = []
        for a in crds:
            if a[2] == prev_a[2]:
                res.append(a)
            else:
                nucs.append(res)
                res = []
            prev_a = a
        return nucs

    # detect strands and store them in separate lists
    def parse_strand(self):
        prev_line = '    1   '
        for line in self.atoms:
            if line[5] != prev_line[5] and line[5] == '1':
                i = self.atoms.index(line)
                self.strand1 = self.atoms[0:i]
                self.strand2 = self.atoms[i:]
            prev_line = line
        return [self.strand1, self.strand2]

    def pair_crds(self):
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

    def printf_pairs(self):
        self.pair_crds()
        with open('pairs_out.dat', 'w') as out:
            for pair in self.base_pairs:
                nuc_a = pair[0]
                nuc_b = pair[1]
                for (a, b) in zip(nuc_a, nuc_b):
                    out.write(('{} {} {} {} {} {}    {} {} {} {} {} {}\n').format(a[0], a[1], a[2], a[3], a[4], a[5], b[0], b[1], b[2], b[3], b[4], b[5]))
                out.write('\n')

    def __init__(self, inp_pdb):
        self.get_atoms(inp_pdb)
        self.mask = ('C3\'', 'O4\'', 'C5\'', 'O3\'', 'P', '05\'')
        self.atoms = self.mask_atoms()
        self.parse_strand()
        self.pair_crds()

inp_pdb = './pdb/3MVA.pdb'
a = dna(inp_pdb)
a.printf_pairs()
