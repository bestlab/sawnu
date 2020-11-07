#!/usr/bin/env python3

# Reference:
# JCP v148, 123329 (2018)
# https://doi.org/10.1063/1.5006954

import sys,numpy,math
from scipy import optimize
from scipy import special

Usage = """

saw_nu.py E N R0

where: E is the mean ratiometric FRET efficiency and N the number of residues,
and R0 is the Foerster radius (in nm!)

"""

def chisq_Egauss(parms, E, R0, Nres):
    R = parms[0]
    rmax = 0.38*Nres
    Nr = 1000
    dr = rmax / float(Nr)
    norm = 0.
    Ecalc = 0.
    for i in range(Nr):
        r=(float(i)+.5)*dr
        Pr = 4*math.pi*r**2*(1.5/(math.pi*R**2))**1.5 * math.exp(-1.5*r**2/R**2)
        norm+= Pr*dr
        Ecalc += Pr*dr/(1+(r/R0)**6) 
    #print(Ecalc, norm)
    return (E-Ecalc)**2

def chisq_Esawnu(parms, E, R0, Nres, b):
    nu = parms[0]
    R = b*Nres**nu
    rmax = 0.38*Nres
    Nr = 1000
    dr = rmax / float(Nr)
    norm = 0.
    Ecalc = 0.
    g=(1.1615-1.)/nu
    delta=1./(1.-nu)
    alpha = (special.gamma((3.+g)/delta)/special.gamma((5.+g)/delta))**(-delta/2.) # hopefully!
    A = 1./( 4*math.pi/(delta)*alpha**(-(3.+g)/delta) * special.gamma((3.+g)/delta) )
    for i in range(Nr):
        r=(float(i)+.5)*dr
        #Pr = 4*math.pi*r**2*(1.5/(math.pi*R**2))**1.5 * math.exp(-1.5*r**2/R**2)
        Pr = A*4.*math.pi/R*(r/R)**(2.+g) * math.exp(-alpha*(r/R)**delta )
        norm+= Pr*dr
        Ecalc += Pr*dr/(1.+(r/R0)**6) 
    #print(Ecalc, norm)
    return (E-Ecalc)**2

def chisq_Esaw(parms, E, R0, Nres):
    nu = 0.588
    R = parms[0]
    rmax = 0.38*Nres
    Nr = 1000
    dr = rmax / float(Nr)
    norm = 0.
    Ecalc = 0.
    g=(1.1615-1.)/nu
    delta=1./(1.-nu)
    alpha = (special.gamma((3.+g)/delta)/special.gamma((5.+g)/delta))**(-delta/2.) # hopefully!
    A = 1./( 4*math.pi/(delta)*alpha**(-(3.+g)/delta) * special.gamma((3.+g)/delta) )
    for i in range(Nr):
        r=(float(i)+.5)*dr
        #Pr = 4*math.pi*r**2*(1.5/(math.pi*R**2))**1.5 * math.exp(-1.5*r**2/R**2)
        Pr = A*4.*math.pi/R*(r/R)**(2.+g) * math.exp(-alpha*(r/R)**delta )
        norm+= Pr*dr
        Ecalc += Pr*dr/(1.+(r/R0)**6) 
    #print(Ecalc, norm)
    return (E-Ecalc)**2


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

if len(sys.argv) != 4:
    sys.stdout.write( Usage )
    sys.exit(0)

EFRET = float(sys.argv[1])
Nres = int(sys.argv[2])
R0_nm = float(sys.argv[3])

init_parm = [R0_nm,]
opt_parm = optimize.fmin(chisq_Egauss, init_parm, (EFRET,R0_nm,Nres))

R_gauss = opt_parm[0]
Rg_gauss = R_gauss/6.**.5

init_parm = [R0_nm,]
opt_parm = optimize.fmin(chisq_Esaw, init_parm, (EFRET,R0_nm,Nres))
R_saw = opt_parm[0]
Rg_saw = R_saw/(6.26)**.5

init_parm = [0.5,]
opt_parm = optimize.fmin(chisq_Esawnu, init_parm, (EFRET,R0_nm,Nres,b_nm))
nu = opt_parm[0]
R_sawnu = b_nm*Nres**nu
gamma = 1.1615
lambda_sawnu = 2*(gamma+2.*nu)*(gamma+2.*nu+1.)/(gamma*(gamma+1))
Rg_sawnu = R_sawnu/lambda_sawnu**.5

sys.stdout.write("================================\n")
sys.stdout.write(" Gaussian chain R = %8.3f\n"%(R_gauss))
sys.stdout.write("Gaussian chain Rg = %8.3f\n"%(Rg_gauss))
sys.stdout.write("            SAW R = %8.3f\n"%(R_saw))
sys.stdout.write("           SAW Rg = %8.3f\n"%(Rg_saw))
sys.stdout.write("         SAW-nu R = %8.3f\n"%(R_sawnu))
sys.stdout.write("        SAW-nu Rg = %8.3f\n"%(Rg_sawnu))
sys.stdout.write("        SAW-nu nu = %8.3f\n"%(nu))
sys.stdout.write("================================\n")



