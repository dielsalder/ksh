#!/bin/bash
#output format:
#frame, phi, theta, residual, radius, pitch, helical_sweep
#the right answer is something like
#1  179.0000    2.0000     0.979477     8.523428    33.493840   726.777015
#

inpfile=B.test.crd
resolution=180
num_atoms=21
num_frames=1
phis=0
phie=180
thes=0
thee=180

./eli55a.o $inpfile $resolution $num_atoms $num_frames $phis $phie $thes $thee
