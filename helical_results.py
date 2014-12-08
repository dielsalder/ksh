    def helical_results(self):
        """calculate helical parameters"""
	x0, y0 = self.center

        thetahelix = 0
        # calculate angle of helix
	for i, crds in enumerate(self.crds):
	    x, y, z = crds
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

        self.zpitch = abs(z_max - z_min)
        self.k = zetapitch * 180/PI / thetahelix
        self.pitch = abs(self.k * 2 * PI)
        return self.pitch
