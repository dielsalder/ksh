from re import search
from collections import Iterable

def flatten(items, ignore_types=(str, bytes, dict)):
    """Flatten nested lists of arbitrary length"""
    for x in items:
        if isinstance(x, collections.Iterable) and not isinstance(x, ignore_types):
            for i in flatten(x):
                yield i
        else:
            yield x

def atom_from_line(line):
    """Read one line of pdb into dict"""
    line = line.split()
    atom = dict(
            atom_num = int(line[1]),
            atom_name = line[2],
            res_name = line[3],
            res_num = int(line[4]),
            x = float(line[5]),
            y = float(line[6]),
            z = float(line[7])
    )
    return atom

def atoms_from_pdb(pdb_name):
    """Read all atoms in pdb file into dict list"""
    with open(pdb_name) as pdb:
        atoms = [atom_from_line(line)
                for line in pdb if line.startswith('ATOM')]
    return atoms

def nucleic_from_pdb(pdb_name):
    """Read only nucleic acid atoms in pdb file into dict list"""
    with open(pdb_name) as pdb:
        atoms = [atom_from_line(line) for line in pdb if re.search
                ('(^ATOM)\s*\S*\s*\S*\s*'
                 '(DA5|DA3|DA|DT5|DT3|DT|DG5|DG3|DG|DC5|DC3|DC)', line)]
    return atoms

def crds_of(atoms):
    """Get xyz coordinates of a subset of atoms"""
    crds = [[atom['x'], atom['y'], atom['z']] for atom in flatten(atoms)]
    return crds

class Dna:
    """Store data for the atoms in one DNA molecule"""
    def __init__(self, pdb_name, mask=[]):
        self.atoms = atoms_from_pdb(pdb_name)
        self.mask = mask
        self.mask_atoms()
        self.group_strands()
        self.base_pairs()
        self.all_crds = crds_of(self.atoms)

    def mask_atoms(self):
        # ugly, fix later
        """Process only atoms specified in mask"""
        if self.mask == []:
            pass
        else:
            masked = [atom for atom in self.atoms if atom['atom_name'] in self.mask]
            self.all_atoms = self.atoms
            self.atoms = masked

    def group_strands(self, n = 1):
        """
        Group atoms by helix strand
        Set n as starf of residue numbering
        """

        prev_atom = {'res_num': n}
        for i, atom in enumerate(self.atoms):
            if (atom['res_num'] != prev_atom['res_num'] and
                    atom['res_num'] == n):
                strand1 = self.atoms[0:i]
                strand2 = self.atoms[i:]
            prev_atom = atom
        for atom in self.atoms:
            if atom in strand1:
                atom['strand'] = 1
            elif atom in strand2:
                atom['strand'] = 2
        return self.atoms

    def base_pairs(self):
        str1 = [atom for atom in self.atoms if atom['strand'] == 1]
        str2 = [atom for atom in self.atoms if atom['strand'] == 2]
        for i, (a, b) in enumerate(zip(str1, str2)):
            a['pair'] = i
            b['pair'] = i

    def get_pairs(self, start = 0, end = 10):
        """Get base pairs"""
        # Both strands are numbered 5' to 3'
        for i in range(start, end):
            yield [atom for atom in self.atoms if atom['pair'] == i]

    def pair_crds(self):
        """Coordinates by base pair"""
        for pair in self.get_pairs():
            yield crds_of(pair)

#    def pdb(self, mask=['P'], out='dna_out.pdb')
#        pass
#        """Write selected atoms to pdb"""
#        include = [atom in flatten(self.atoms) if atom['atom_name'] in mask]
