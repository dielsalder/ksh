import re
import operator
import itertools
import collections

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

class DnaMolecule:
    """Store data for the atoms in one DNA molecule"""
    def __init__(self, pdb_name, mask=[]):
        self.atoms = atoms_from_pdb(pdb_name)
        self.mask = mask
        self.mask_atoms()
        self.group_strands()
        self.group_nucs()
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

    def group_strands(self):
        """Group atoms by helix strand"""
        n = 1   # Residue numbering start
        prev_atom = {'res_num': n}
        for i, atom in enumerate(self.atoms):
            if (atom['res_num'] != prev_atom['res_num'] and
                    atom['res_num'] == n):
                strand1 = self.atoms[0:i]
                strand2 = self.atoms[i:]
                self.atoms = [strand1, strand2]
            prev_atom = atom
        #for strand in self.atoms:
        #    # sort by residue
        #    strand.sort(key=operator.itemgetter('res_num'))
        return self.atoms

    def group_nucs(self):
        """Group atoms within strand by nucleotide"""
        keyfunc = operator.itemgetter('res_num')
        self.atoms = [[list(grp) for key, grp in itertools.groupby(sorted(strand, key=keyfunc),
            key=keyfunc)] for strand in self.atoms]
        return self.atoms

    def pair_nucs(self, n=0):
        """Get base pairs"""
        # Both strands are numbered 5' to 3'
        for (a, b) in zip(self.atoms[0], reversed(self.atoms[1])):
            yield (a, b)

    def pair_crds(self):
        """Coordinates by base pair"""
        for pair in self.pair_nucs():
            yield crds_of(pair)
