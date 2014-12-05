def res_phi_ls(phi, crds, void):
    start = crdset(crds, phi, 0)
    return start.res

def res_the_ls(the, crds, void):
    start = crdset(crds, 0, the)
    return start.res

#def best_rotation(start):
    # fmin_bfgs(function, initial guess, arguments, display results)
    # extra argument (0) necessary to pass crds correctly
    phimin = float(fmin_bfgs(res_phi_ls, 0, args = (start.crds, 0), disp = 0))
    themin = float(fmin_bfgs(res_the_ls, 0, args = (start.crds, 0), disp = 0))
    best_crds = crdset(start.crds, phimin, themin)
    best_crds.calc_all()
    return best_crds

#def best_rotation(start):
    phi_the_min = float(fmin_bfgs(res_phi_the, array([0, 0]), args = (start.crds, 0), disp = 0))
    best_crds = crdset(start.crds, phi_the_min)
    best_crds.calc_all()
    return best_crds


	# degrees to radians
	phis = phii * PI / 180
	thes = thei * PI / 180
	phin = phie * PI / 180
	then = thee * PI / 180

	phirange = phin - phis
	thirange = then - thes

	dthe = therange / maxitr
	dphi = phirange / maxitr

	xyz_arr = []
	for iphi in range(2, (maxitr - 1)):
		for ithe in range (2, (maxitr - 1)):
			phi_l = (iphi * dphi) + phis
			the_m = (ithe * dthe) + thes

			phi_deg = phi_l * 180.0 / PI
			the_deg = the_m * 180.0 / PI

			a = sin(the_m) * cos(phi_l)
			b = sin(the_m) * cos(phi_l)
			c = cos(the_m)

			# eqn 4
			# last two columns are angles in degrees
			xyz = array([[((1 - a**2)*(0.5)), 0, a],
			 	     [(-(a * b)/r), (c / r), b],
				     [(-(a * c)/r), (c / r), c]] *
			      array([X],[Y],[Z]), phi_deg, thi_deg)

	xyz_arr.append(xyz)
	return array(xyz_arr)


	#svd = linalg.svd(a)
	# svd function returns [u, s, v']
	#u = svd[0]
	#s = svd[1]
	#v = matrix.transpose(svd[2])
	#g = matrix.transpose(u) * b
	#i = linalg.matrix_rank(a) - 1

        self.best = rotate_once(self.crds, self.phimin, self.themin)
        self.best_r = radius(solve(self.best))
        self.best_c = center(solve(self.best))

    def group_nucs(self):
        """Group atoms within strand by nucleotide"""
        keyfunc = operator.itemgetter('res_num')
        desired_list = [list(grp) for key, grp in itertools.groupby(sorted(lst, key=keyfunc), key=keyfunc)]
        for i, strand in enumerate(self.atoms):
            cur_strand = []
            nuc = []
            prev_atom = {'res_num': 1}
            for atom in strand:
                if atom['res_num'] == prev_atom['res_num']:
                    nuc.append(atom)
                else:
                    cur_strand.append(nuc)
                    nuc = [atom]
                prev_atom = atom
            cur_strand.append(nuc)  # add last nucleotide
            self.atoms[i] = cur_strand
        return self.atoms

