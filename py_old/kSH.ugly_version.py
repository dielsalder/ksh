from numpy import *
from scipy.optimize import *

PI = 3.14159265359

inpcrds = 'bdna_rotated.dat'
with open(inpcrds) as f:
    crds = loadtxt(f)

phis_deg = 0.1
phin_deg = 180
thes_deg = 0.1
then_deg = 180
maxitr = 50

# convert degrees to radians
phis = radians(phis_deg)
phin = radians(phin_deg)
thes = radians(thes_deg)
then = radians(then_deg)

# get next angle in range
def new_ang(angs, angn, maxitr):        # start, end, maxitr
        angr = angs
        angstep = (angn - angs) / maxitr
        while angr <= angn:
                yield angr
                angr += angstep

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
        return svd[0]   # x

def radius(x):
        return (x[1]**2 + x[2]**2 + x[0])**(0.5)

def center(x):
        return [x[1], x[2]]

def pitch(crds, x):
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

# calculate residual using x
def res_x(crds, x):
        res = 0
        x0 = x[1]
        y0 = x[2]
        r = radius(x)
        for i in crds:
                x = i[0]
                y = i[1]
                res += abs((x - x0)**2 + (y - y0)**2 +
                            r**2 - (2 * r) *
                           ((x - x0)**2 + (y - y0) ** 2)**0.5)
        return res

# residual as function of phi and theta
# rotates and then solves res_x using rotated coordinates
def res_phi_the(phi, the, old_crds):
        r_crds = rotate_once(old_crds, phi, the)
        x = solve(r_crds)
        return res_x(r_crds, x)

# residual as function of phi
def res_phi(phi, old_crds, the):
        r_crds = rotate_once(old_crds, phi, the)
        res = res_x(r_crds, solve(r_crds))
        return res

# residual as function of theta
def res_the(the, old_crds, phi):
        r_crds = rotate_once(old_crds, phi, the)
        res = res_x(r_crds, solve(r_crds))
        return res

def res_lstsq(crds):
        return abs(solve(crds)[1])

def res_phi_ls(phi, old_crds, the):
        r_crds = rotate_once(old_crds, phi, the)
        res = res_lstsq(r_crds)
        return res

def res_the_ls(the, old_crds, phi):
        r_crds = rotate_once(old_crds, phi, the)
        res = res_lstsq(r_crds)
        return res

class crdset:
        def __init__(self, crds):
                self.x = solve(crds)
                self.res = res_x(crds, self.x)
                self.radius = radius(self.x)
                self.center = center(self.x)
                self.phimin = float(fmin_bfgs(res_phi_ls, 0, args = (crds, 0), disp = 0))
                self.themin = float(fmin_bfgs(res_the_ls, 0, args = (crds, self.phimin), disp = 0))
                self.best = rotate_once(crds, self.phimin, self.themin)
                self.best_r = radius(solve(self.best))
                self.best_c = center(solve(self.best))

inpcrdset = crdset(crds)
print inpcrdset.phimin, inpcrdset.themin
print inpcrdset.best_r, inpcrdset.best_c
# unrotated coordinates
print "\nOriginal coordinates"
print "Radius: ", radius(solve(crds)), "\tPitch: ", pitch(crds, solve(crds))
print "Residual: ", res_x(crds, solve(crds)),

print "\nResidual: ", res_x(rotate_once(crds, 0, 0), solve(rotate_once(crds, 0, 0)))


# minimize res(phi), then res(the)
phimin = float(fmin_bfgs(res_phi, 0, args = (crds, 0), disp = 0))
themin = float(fmin_bfgs(res_the, 0, args = (crds, phimin), disp = 0))

x_min = solve(rotate_once(crds, phimin, themin))
res_min = res_phi_the(phimin, themin, crds)

# optimized circle facts
r_min = radius(x_min)
c_min = center(x_min)
pitch_min = pitch(crds, x_min)

print "\nMinimizer"
print "x = ", x_min
print "phi = ", phimin, "\tthe = ", themin
print "Radius: ", r_min, "\tCenter: ", c_min, "\nResidual: ", res_min, #"\nPitch: ", pitch_min

# determine lowest res by minloc in array
#residuals = []
#phithe = []
#for phi in new_ang(phis, phin, maxitr):
#        for the in new_ang(thes, then, maxitr):
#                new_crds = rotate_once(crds, phi, the)
#                x = solve(new_crds)
#                residuals.append(res_x(new_crds, x))
#                phithe.append([phi, the])
#resmin = min(residuals)
#phithemin = phithe[residuals.index(resmin)]
#phimin = phithemin[0]
#themin = phithemin[1]
#xmin = solve(rotate_once(crds, phimin, themin))
#print "\nMinimum in array"
#print "phi = ", phimin, "\tthe = ", themin, "\nResidual: ", resmin
#print "Radius: ", radius(xmin)
#
