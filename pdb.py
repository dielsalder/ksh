import collections
import ksh

def flatten(x, ignore_types=(str, bytes, dict)):
    for i in x:
        if isinstance(x, collections.Iterable) and not isinstance(x, ignore_types):
            for i in flatten(x):
                yield i
        else:
            yield x

def atom(line):
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

def get_crds(atoms):
    """Get xyz coordinates of selected atoms"""
    crds = [[a['x'], a['y'], a['z']] for a in atoms]
    return crds

class Pdb:
    """pdb handler"""
    def __init__(self, name):
        self.filename = name

    def get_all(self):
        """Read all atoms in pdb file"""
        with open(self.filename) as pdb:
            atoms = [atom(line)
                    for line in pdb if line.startswith('ATOM')]
        return atoms

    def get_nucleic(self):
        """Read only nucleic acid atoms in pdb file"""
        with open(self.filename) as pdb:
            atoms = [atom(line) for line in pdb if search
                    ('(^ATOM)\s*\S*\s*\S*\s*'
                     '(DA5|DA3|DA|DT5|DT3|DT|DG5|DG3|DG|DC5|DC3|DC)', line)]
        return atoms

class Molecule:
    """pdb file"""
    def __init__(self, atoms):
        self.atoms = atoms
        self.fits = []

    def select(self, res = [], atom_name = []):
        """Select atoms by residue and element"""
        selection = (
            [atom for atom in self.atoms if atom['res_name'] in res
                and atom['atom_name'] in atom_name])
        return selection

    def get_res(self, res_num):
        return [a for a in self.atoms if a['res_num'] == res_num]

    def write(self, atoms, out = open('atoms.pdb', 'w')):
        """Write selected atoms to pdb"""
        out.write('REMARK      generated by pdb.py\n')
        for atom in atoms:
            vals = (['ATOM', atom['atom_num'], atom['atom_name'],
                atom['res_name'], atom['res_num'],
                atom['x'], atom['y'], atom['z'],
                '1.00', '0.00', '\n'])
            line = '    '.join(str(v) for v in vals)
            out.write(line)

    def minimize(self, selection):
        """fit minimized helix"""
        crds = get_crds(selection)
        fit = Fit(crds)
        min_fit = (fit.minimize())
        self.fits.append(min_fit)
        return min_fit

class Dna(Molecule):
    """DNA molecule (from a pdb file)"""
    def __init__(self, atoms):
        self.atoms = atoms
        self.strands()
        self.pairs()
        self.all_crds = get_crds(self.atoms)

    def strands(self, n = 1):
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

    def get_strand(self, strand):
        return [atom for atom in self.atoms if atom['strand'] == strand]

    def pairs(self):
        str1 = self.get_strand(1)
        str2 = self.get_strand(2)
        current = 0
        for (a, b) in zip(str1, str2):
            if not a['res_num'] == current:
                current += 1
            a['pair'] = current
            b['pair'] = current
        # kludge: don't know why, but with test file ideal_bdna.pdb atoms
        # in strand 1 starting at 345 are unpaired
        # this assigns last pair number to all unpaired atoms
        for a in self.atoms:
            if not 'pair' in a:
                a['pair'] = current
        return self.atoms

    def get_pair(self, i_pair):
        """Get one base pair"""
        return [a for a in self.atoms if a['pair'] == i_pair]

    def get_pairs(self, start, end):
        """Get base pairs from start to end, inclusive"""
        # Both strands are numbered 5' to 3'
        return [a for a in self.atoms if a['pair'] in range(start, end + 1)]

class Fit:
    """Results of one fit"""
    def __init__(self, atoms):
        self.crds = get_crds(atoms)
        self.crdset = ksh.crdset(self.crds)

    def no_rotate(self):
        """Fit without rotation"""
        fit = self.crdset.calc_all()
        self.no_rotate = fit
        return self.no_rotate

    def minimize(self):
        """Minimize rotation using ksh's best_rotation"""
        rotation = ksh.best_rotation(self.crdset)
        best = rotation.calc_all()
        self.res = best[1]
        self.best = best
        return self.best

    def write_minimize(self, Molecule):
        """Write minimized fit to Molecule.fits"""
        Molecule.fits.append(self.best)
        return Molecule.fits
