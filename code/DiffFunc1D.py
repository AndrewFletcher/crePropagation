# Functions used for cosmic ray diffusion code
# periodic, 1D Cartesian coordinates

import numpy as np

def courantDiffusion(D, h, dt):
    # check stability condition
	cond = 2*D*dt/(h**2)
	if cond < 1.0:
	  print "Courant Diffusion: ", str(cond), "PASS"
	if cond >= 1.0:
	  print "Courant Diffusion: ", str(cond), "FAIL"

def initN(rmin, rmax, nr, n0, distribution):
    # source spatial distribution
    n = np.linspace(rmin, rmax, nr)
    if distribution == "triangle":
        midPoint = nr/2.0
        for i in range(nr):
            if i <= midPoint:
                n[i] = (2*n0/(rmax-rmin))*(n[i]-rmin)
            else:
                n[i] = (2*n0/(rmax-rmin))*(n[nr-1]-n[i])
        return n
    if distribution == "delta":
        n[:] = 0
        n[rmin] = n0
        return n
    if distribution == "twoArm":
        arm1Centre = (rmax-rmin)/4
        arm2Centre = 3*(rmax-rmin)/4
        width = (rmax-rmin)/20
        n=(n0*np.exp(-(n-arm1Centre)**2/(2*width**2)) 
            + n0*np.exp(-(n-arm2Centre)**2/(2*width**2)))
        return n

def diff(n, nr, dr):
    # periodic first derivative
	dn = (np.roll(n,-1) - np.roll(n,1))/dr
	return dn

def diff2(n, nr, dr):
    # periodic second derivative
	d2n = (np.roll(n,-1) - 2*n + np.roll(n,1))/(dr**2)
	return d2n

def rhsDiff(n, nr, dr):
    # RHS of diffusion term
    rhs = diff2(n, nr, dr)
    return rhs

def courantEnergy(Phi, dt, deltaE):
    # check stability condition
    stabE = Phi*dt/deltaE[0]
    print "Energy Courant = ", stabE
    if stabE>1.0:
        print "Courant condition not satisfied"
        raise SystemExit

def initE(E, n0, gamma):
    # source power-law distribution in E
    return n0*E**(gamma)
    
def Ebc(N,E):
    # find N at ghost zone energy
    gamma1=np.log(N[0]/N[1])/np.log(E[1]/E[2])  # E[0] is ghost energy
    N0=np.exp(np.log(N[0])+gamma1*np.log(E[0]/E[1]))
    return N0

def rhsE(Phi, N, E, deltaE):
    # RHS of energy-loss term
    ghostN = Ebc(N[0:2],E[0:3])  # need N0 to go with E0
    tempN = np.append(ghostN,N)  # original N starts at E1
    diff = (tempN-np.roll(tempN,1))  # 1st order upwind difference
    upN = diff[1:]/deltaE  # upwind derivative, tempN[0] is ghost
    rhs=Phi*(upN*E[1:]**2 + 2*E[1:]*N)
    return rhs

def floorN(N):
    # put a floor on N: i.e. if N<0, N=nan
    N[N<0.0]=np.nan
    return N

def specInd(N1,N2,E1,E2):
    nN = np.shape(N2)[0]
    gamma=np.zeros(nN)
    for i in range(nN):
        if N2[i]>0:
            gamma[i] = -1.0*(np.log(N1[i]/N2[i]) / np.log(E1/E2))
        else:
            gamma[i]=np.nan
    return gamma

