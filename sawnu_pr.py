#!/usr/bin/env python3

# Reference:
# JCP v148, 123329 (2018)
# https://doi.org/10.1063/1.5006954

import sys,numpy,math
from scipy import optimize
from scipy import special

Usage = """

saw_nu_pr.py nu N

where: nu is the scaling exponent and N the number of residues

"""


def pr_sawnu_find_a_alpha(R,nu):
    # pre-solve for A, alpha so we don't need to compute them each time Pr is calculated
    g=(1.1615-1.)/nu
    delta=1./(1.-nu)
    alpha = (special.gamma((3.+g)/delta)/special.gamma((5.+g)/delta))**(-delta/2.) # hopefully!
    A = 1./( 4*math.pi/(delta)*alpha**(-(3.+g)/delta) * special.gamma((3.+g)/delta) )
    return A, alpha

def pr_sawnu(r,R,nu,A,alpha):
    g=(1.1615-1.)/nu
    delta=1./(1.-nu)
    Pr = A*4.*math.pi/R*(r/R)**(2.+g) * math.exp(-alpha*(r/R)**delta )
    return Pr

b_nm = 0.55 # fixed b for now

if len(sys.argv) != 3:
    sys.stdout.write( Usage )
    sys.exit(0)

nu = float(sys.argv[1])
Nres = int(sys.argv[2])

R = b_nm * Nres**nu

nstep = 500
max_r = R*10.
dr = max_r/float(nstep)

A,alpha = pr_sawnu_find_a_alpha(R,nu)

S=0.
for s in range(nstep):
    r = dr*float(s)
    pr = pr_sawnu(r,R,nu,A,alpha)
    S+=pr*dr
    sys.stdout.write("%12.6f %12.6e %12.6e\n"%(r,pr,S))



