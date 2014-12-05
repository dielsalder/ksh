from numpy import *
from scipy.optimize import *

PI = 3.14159265359

inpcrds = 'bdna.pdb'
with open(inpcrds) as inp:
    crds = loadtxt(inp)

phi_start_deg = 0.1
phi_end_deg = 180
the_start_deg = 0.1
the_end_deg = 180
maxitr = 50

# convert degrees to radians
phi_start = radians(phi_start_deg)
phi_end = radians(phin_end_deg)
the_start = radians(the_start_deg)
the_end = radians(the_end_deg)

# get next angle in range
def new_ang(start, end, maxitr):
    new = start
    step = (end - start) / maxitr
    while new <= end:
        yield new
        new += step

# rotate xyz coordinates by angles phi and theta
def rotate_once(crds, phi, the):
    new_crds = []
    a = sin(the) * cos(phi)
    b = sin(the) * cos(phi)
    c = cos(the)
    r = (1-(a**2))**(0.5)

    for iatom in crds:
        new_xyz = dot(array([[((1 - a**2)*(0.5)), 0, a],
                            [(-(a * b)/r),  (c / r), b],
                            [(-(a * c)/r), -(b / r), c]]),
                          iatom)
        new_xyz[0] = new_xyz[0] * 2
        new_crds.append(new_xyz)
    return array(new_crds)

# build a matrix
def build_amat(crds):
    amat = insert(crds, 0, values = 1, axis = 1)
    return amat[:, [0, 1, 2]]

# build b matrix
def build_bmat(crds):
    bmat = []
    for iatom in crds:
        bmat.append(iatom[0]**2 + iatom[1]**2)
    return array(bmat)

# solve ax = b
def solve(crds):
    a = build_amat(crds)
    b = build_bmat(crds)
    svd = linalg.lstsq(a, b)
    # linalg.lstsq returns
    #       [x, sum of residuals, rank(a), singular values of a]
    return svd
    #return svd[0]   # x

# residual according to least squares routine
def res_lstsq(crds):
    return solve(crds)[1]

def pitch(crdset):
    crds = crdset.crds
    x = crdset.x

    x0 = x[1]
    y0 = x[2]

    thetahelix = 0
    # calculate angle of helix
    for i in range(1, len(crds)-1):
        x = crds[i, 0]
        y = crds[i, 1]
        vecix = x - x0  # ith x
        veciy = y - y0  # ith y
        vecx1 = crds[i + 1, 0] - x0     #i+1th x
        vecy1 = crds[i + 1, 1] - x0     #i+1th y

        doti = vecix * vecx1 + veciy * vecy1
        magi = (vecix**2 + veciy**2)**0.5       # magnitudes
        magi1 = (vecx1**2 + vecy1**2)**0.5

        dotovermags = doti / (magi * magi1)
        thetastep = arccos(dotovermags)
        thedstep = thetastep * 180/PI

        thetahelix += thedstep

    # calculate kappa
    z_max = max(crds[:,2])
    z_min = min(crds[:,2])

    zetapitch = abs(z_max - z_min)
    kappa = zetapitch * 180/PI / thetahelix
    pitch = abs(kappa * 2 * PI)
    return pitch

def res_phi_ls(phi, crdset):
    crdset.rotate(phi, 0)
    return crdset.res

def res_the_ls(the, crdset);
    crdset.rotate(0, the)
    return crdset.res

def find_best_rotation(old_crds):
    phimin = float(fmin_bfgs(res_phi_ls, 0, args = (self.old_crds, 0), disp = 0))
    themin = float(fmin_bfgs(res_the_ls, 0, args = (self.old_crds, phimin), disp = 0))
    best_crds = crdset(old_crds.crds, phimin, themin)
    return best_crds

class crdset:
    def lstsq(self):
        lstsq = solve(self.crds)
        self.x = lstsq[0]
        self.res = lstsq[1]
        x = self.x
        self.radius = x[1]**2 + x[2]**2 + x[0])**(0.5)
        self.center = [x[1], x[2]]
    def rotate(self, phi, the):
        self.crds = rotate_once(self.crds, phi, the)
        self.lstsq()
    def __init__(self, crds, *phi_the):
        self.crds = crds
        if phi_the:
            self.phi = phi_the[0]
            self.the = phi_the[1]
            self.rotate(crds, phi, the)
        else:
            self.lstsq()

# print results of calculation with optional header
def print_out(crdset, *args):
    print "\n", args
    print "\nresidual: ", crdset.res
    print "\nRadius: ", crdset.radius
    print "\nCenter: ", crdset.center
    print "\nHelical pitch: ", crdset.pitch

print_out(crds, "Original coordinates")
best = find_best_rotation(crds)
print_out(best, "Best fit")
