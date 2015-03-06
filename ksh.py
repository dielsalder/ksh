#from numpy import *
# if name not defined, it's in numpy - continue tomorrow
import numpy as np
from numpy import sin, cos, tan, dot, arccos
from scipy.optimize import fmin_bfgs

PI = 3.14159265359
# numpy ignores divide by zero error
# doesn't work but this is apparently how it should work?
np.seterr(divide = 'ignore')

# linear least squares problem for circle-fitting
def build_amat(crds):
    """build a matrix"""
    amat = np.insert(crds, 0, values = 1, axis = 1)
    return amat[:, [0, 1, 2]]

def build_bmat(crds):
    """build b matrix"""
    bmat = []
    for iatom in crds:
        bmat.append(iatom[0]**2 + iatom[1]**2)
    return np.array(bmat)

def solve(crds):
    """solve ax = b"""
    a = build_amat(crds)
    b = build_bmat(crds)
    svd = np.linalg.lstsq(a, b)
    # linalg.lstsq returns
    #       [x, sum of residuals, rank(a), singular values of a]
    return svd
    #return svd[0]   # x

def res_phi_the(phi_the, crds, void):
    """residual as function of phi and theta"""
    phi, the = phi_the
    start = crdset(crds, phi, the)
    return start.res

def resv_phi_the(phi_the, crds, void):
    """other residual as function of phi and theta"""
    phi, the = phi_the
    start = crdset(crds, phi, the)
    return start.res_v

class crdset:
    def lstsq(self):
        """Solve ax = b"""
        lstsq = solve(self.crds)
        self.x = lstsq[0]
        self.calc_res_v()
        # get residual from least squares solution
        self.res_svd = float(lstsq[1])
        self.singvals = lstsq[3]
        return self.lstsq

    def lstsq_results(self):
        """Circle-fitting data"""
        x = [i + 0j for i in self.x]
        self.radius = x[1]**2 + x[2]**2 + x[0]**0.5
        self.center = [x[1], x[2]]
        return [self.radius, self.center]

    def calc_res_v(self):
        """Vageli's residual"""
        x = self.x
        res = 0
        x0 = x[1]
        y0 = x[2]
        self.lstsq_results()
        r = self.radius
        for i in self.crds:
                x = i[0]
                y = i[1]
                res += abs((x - x0)**2 + (y - y0)**2 +
                            r**2 - (2 * r) *
                           ((x - x0)**2 + (y - y0) ** 2)**0.5)
        self.res_v = res
        return self.res_v

    def helical_results(self):
        """calculate helical parameters"""
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
            thedstep = thetastep * 180/PI

            thetahelix += thedstep

        # calculate kappa
        # max and min of z-coordinate column
        z_max = max(self.crds[:,2])
        z_min = min(self.crds[:,2])

        zetapitch = abs(z_max - z_min)    # height of helix
        self.k = zetapitch * 180/PI / thetahelix
        self.pitch = abs(self.k * 2 * PI)
        self.zpitch = zetapitch
        return self.pitch

    def rotate(self, phi, the):
        """rotate by phi and theta"""
        new_crds = []
        a = sin(the) * cos(phi)
        b = sin(the) * cos(phi)
        c = cos(the)
        r = (1-(a**2))**(0.5)

        for iatom in self.crds:
            new_xyz = dot(np.array([[((1 - a**2)*(0.5)), 0, a],
                                [(-(a * b)/r),  (c / r), b],
                                [(-(a * c)/r), -(b / r), c]]),
                              np.array(iatom))
            new_xyz[0] = new_xyz[0] * 2
            new_crds.append(new_xyz)
        self.crds = np.array(new_crds)
        # recalculate residual, necessary for minimization
        self.lstsq()
        return self.crds

    def iter_rotate(self):
        best_res = 1000000000000000000000000000
        for phi in range(0, 360):
            for the in range(0, 360):
                self.rotated_crds = self.rotate(phi, the)
                self.lstsq()
                if self.res_svd < best_res:
                    best_res = self.res_svd
                    best_phi = phi
                    best_the = the
        self.phi = best_phi
        self.the = best_the
        self.best_rotate = (best_phi, best_the)
        self.best_res = best_res
        return [self.rotate_crds, best_phi, best_the]

    def calc_all(self):
        """Calculate everything"""
        #self.lstsq()            # update least-squares stuff
        #self.calc_res_v()
        self.lstsq_results()    # circle parameters
        self.helical_results()  # helical parameters
        # uses svd residual
        # simplify phi, the
        return [self.lstsq(), self.res_svd, self.lstsq_results,
            self.helical_results(), (self.phi, self.the)]

    def print_results(self, title=''):
        """Print results of calculations with optional title"""
        print "\n", title
        try:
            print "phi = ", self.phi, "\tthe = ", self.the
        except AttributeError:
            pass
        print "residual: ", self.res
        print "residual 2: ", self.res_v
        print "x = \t", self.x
        print "Radius: ", self.radius
        print "Center: ", self.center
        print "Singular values:", self.singvals
        print "zetapitch\t", self.zpitch
        print "K = ", self.k
        print "Helical pitch: ", self.pitch

    def __init__(self, crds, *phi_the):
        self.crds = np.array(crds)
        if phi_the:
            self.phi, self.the = phi_the
            self.rotate(self.phi, self.the)
        else:
           self.lstsq()

# use scipy's bfgs optimization function
def best_rotation(start):
    phi_the_min = fmin_bfgs(resv_phi_the, np.array([0, 0]), args = (start.crds, 0), disp = 0)
    phi, the = phi_the_min
    best_crds = crdset(start.crds, phi, the)
    best_crds.calc_all()
# best_rotation should be in crdset class, phi/the annoying to access - fix this
    return (best_crds, phi, the)

def test_fit(reference, rotated):
    reference.calc_all()
    print_results(reference, "True best coordinates")
    rotated.calc_all()
    print_results(rotated, "Best fit")
