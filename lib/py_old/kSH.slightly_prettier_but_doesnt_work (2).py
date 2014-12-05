from numpy import *
from scipy.optimize import *

PI = 3.14159265359

# linear least squares problem for circle-fitting
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

# residual as function of phi and theta
def res_phi_the(phi_the, crds, void):
    phi = phi_the[0]
    the = phi_the[1]
    start = crdset(crds, phi, the)
    return start.res

class crdset:
    def lstsq(self):
        lstsq = solve(self.crds)
        self.x = lstsq[0]
#       self.calc_res_v()
        # get residual from least squares solution
        self.res = float(lstsq[1])
        self.singvals = lstsq[3]
        return self.lstsq

    # calculate circle-fitting data
    def lstsq_results(self):
        x = self.x
        self.radius = x[1]**2 + x[2]**2 + x[0]**0.5
        self.center = [x[1], x[2]]
        return [self.radius, self.center]

    # Vageli residual
    def calc_res_v(self):
        x = self.x
        res = 0
        x0 = x[1]
        y0 = x[2]
        self.lstsq_results()
        r = self.radius
        for i in crds:
                x = i[0]
                y = i[1]
                res += abs((x - x0)**2 + (y - y0)**2 +
                            r**2 - (2 * r) *
                           ((x - x0)**2 + (y - y0) ** 2)**0.5)
        self.res_v = res
        return self.res_v

    # helical parameters
    def helical_results(self):
        x0 = self.center[0]
        y0 = self.center[1]

        thetahelix = 0
        # calculate angle of helix
        for i in range(0, len(self.crds)-1):
            x = self.crds[i, 0]
            y = self.crds[i, 1]
            vecix = x - x0  # ith x
            veciy = y - y0  # ith y
            vecx1 = self.crds[i + 1, 0] - x0     #i+1th x
            vecy1 = self.crds[i + 1, 1] - y0     #i+1th y

            doti = vecix * vecx1 + veciy * vecy1
            magi = (vecix**2 + veciy**2)**0.5       # magnitudes
            magi1 = (vecx1**2 + vecy1**2)**0.5

            dotovermags = doti / (magi * magi1)
            thetastep = arccos(dotovermags)
            thedstep = thetastep * 180/pi

            thetahelix += thedstep

        # calculate kappa
        z_max = max(self.crds[:,2])
        z_min = min(self.crds[:,2])

        zetapitch = abs(z_max - z_min)
        self.zpitch = zetapitch
        kappa = zetapitch * 180/PI / thetahelix
        self.k = kappa
        pitch = abs(kappa * 2 * PI)
        self.pitch = pitch
        return self.pitch

    # rotate by phi and theta
    def rotate(self, phi, the):
        new_crds = []
        a = sin(the) * cos(phi)
        b = sin(the) * cos(phi)
        c = cos(the)
        r = (1-(a**2))**(0.5)

        for iatom in self.crds:
            new_xyz = dot(array([[((1 - a**2)*(0.5)), 0, a],
                                [(-(a * b)/r),  (c / r), b],
                                [(-(a * c)/r), -(b / r), c]]),
                              iatom)
            new_xyz[0] = new_xyz[0] * 2
            new_crds.append(new_xyz)
        self.crds = array(new_crds)
        # recalculate residual, necessary for minimization
        self.lstsq()
        return self.crds

    # calculate everything
    def calc_all(self):
        self.lstsq()            # update least-squares stuff
        self.lstsq_results()    # circle parameters
        self.helical_results()  # helical parameters

    def __init__(self, crds, *phi_the):
        self.crds = array(crds)
        if phi_the:
            self.phi = phi_the[0]
            self.the = phi_the[1]
            self.rotate(self.phi, self.the)
        else:
           self.lstsq()

# use scipy's bfgs optimization function
def best_rotation(start):
    phi_the_min = fmin_bfgs(res_phi_the, array([0, 0]), args = (start.crds, 0), disp = 0)
    phi = phi_the_min[0]
    the = phi_the_min[1]
    best_crds = crdset(start.crds, phi, the)
    best_crds.calc_all()
    return best_crds

#  print results of calculations with header
def print_results(crdset, name, *args):
    print name
    try:
        print "phi = ", crdset.phi, "\tthe = ", crdset.the
    except AttributeError:
        pass
#   print "x = ", crdset.x
    print "residual: ", crdset.res
    print "Radius: ", crdset.radius
    print "Center: ", crdset.center
    print "Singular values:", crdset.singvals
    print "zetapitch\t", crdset.zpitch
    print "K = ", crdset.k
    print "Helical pitch: ", crdset.pitch

def test_fit(reference, rotated):
    reference.calc_all()
    print_results(reference, "True best coordinates")
    rotated.calc_all()
    print_results(rotated, "Best fit")

with open('B.test.crd') as inpcrds:
    crds = loadtxt(inpcrds)
    start_crds = crdset(crds)
with open('bdna_rotated.dat') as inpcrds:
    crds = loadtxt(inpcrds)
    rotated = crdset(crds)
with open('bdna.pdb') as inpcrds:
    crds = loadtxt(inpcrds)
    bdna_crds = crdset(crds)

start_crds.calc_all()
print_results(start_crds, "results")
