import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import DiffFunc2D as df

# params
D=10.0
rmin=0.0; rmax=10.0; nr=101; dr=1.0;
n0 = 10.0; 
Emin=-2; Emax=1; nE = 100; gamma = -2.1; 
Phi = 3.0
nt=100; dt=0.00001; 
escT = 0.003  # escape time in time units (same as dt)
source = "yes"
sourceTime = -1.0  # how many times source will pass through domain
record = "no"

# set up arrays
history = np.zeros((nr, nE-1, nt))  # nE-1 as nE[0] is ghost
r = np.linspace(rmin, rmax, nr)
E = np.logspace(Emin,Emax,nE)  # includes the ghost zone, E_min
deltaE = E[1:]-E[:-1]  # the delta Es are [E2-E1,E3-E2,...]
t = np.linspace(0,nt*dt,nt)  # time
nTotal = np.zeros((nt))  # store total number of CRE
nEmax = np.zeros((nt)) # store N(:,Emax)
nEmid = np.zeros((nt)) #       N(:,Emid)
nEmin = np.zeros((nt)) #       N(:,Emin)

# test stability
df.courantDiffusion(D, dr, dt)
df.courantEnergy(Phi, dt, deltaE)

# source distribution
NN0 = df.initN(rmin, rmax, nr, n0, "twoArm")
# Energy distribution
NE0 = df.initE(E[1:], NN0, gamma)  # excludes the ghost zone N(E_min)
n = NE0
nTotal[0] = np.sum(n)

if record == "yes":
    history[:,:,0] = n

i=1
while i<nt:
    # first the diffusion
    k1 = dt*D*df.rhsDiff(n, nr, dr)
    nNew = n+k1/2.0
    k2 = dt*D*df.rhsDiff(nNew, nr, dr)
    n = n+k2
    # now energy losses
    k1 = dt*df.rhsE(Phi, n, E, deltaE)
    nNew = n+k1/2.0
    k2 = dt*df.rhsE(Phi, nNew, E, deltaE)
    n = n+k2
    
    n = n*(1.0-(dt/escT))  # escape losses
 
    if source == "yes":
        if sourceTime != 0.0:
            step = int(np.floor(i/(nt/(sourceTime*nr))))
        if sourceTime == 0.0:
            step=0
        n = n + np.roll(NE0,step, axis=0)

    nTotal[i] = np.sum(n)

    if record == "yes":
        history[:,:,i] = n
   
    i = i+1

# calculate the spectral index
specLoMid = df.specInd(n[:,0],n[:,nE/2],E[1],E[1+nE/2])
specLoHi = df.specInd(n[:,0],n[:,nE-2],E[1],E[nE-1])

# plotting params
fs=20
lw=2

fig = plt.figure(figsize=(8,10))
fig.subplots_adjust(hspace=0.8, wspace=0.5)

ax1 = fig.add_subplot(311)
#ax1.set_title(r"$N(x,[Emax,Emid,Emin])$")
ax1.set_xlabel(r"$\phi \ [2\pi]$",fontsize=fs)
ax1.set_ylabel(r"$N$",fontsize=fs)
ax1.set_yscale('log')
ax1.plot(r/rmax,n[:,0],'b',linewidth=lw)
ax1.plot(r/rmax,n[:,nE/2],'b-.',linewidth=lw)
ax1.plot(r/rmax,n[:,nE-2],'b:',linewidth=lw)

ax2 = fig.add_subplot(323)
ax2.set_title("Arm")
ax2.set_xlabel(r"$E$",fontsize=fs)
ax2.set_ylabel(r"$N$",fontsize=fs)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.plot(E[1:],n[nr/4,:],'b',linewidth=lw)

ax3 = fig.add_subplot(324)
ax3.set_title("Interarm")
ax3.set_xlabel(r"$E$",fontsize=fs)
ax3.set_ylabel(r"$N$",fontsize=fs)
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.plot(E[1:],n[nr/2,:],'b',linewidth=lw)

#~ ax4 = fig.add_subplot(313)
#~ ax4.set_title("Total N")
#~ ax4.set_xlabel(r"$t$")
#~ ax4.set_ylabel(r"$N$")
#~ ax4.plot(t,nTotal)

ax4 = fig.add_subplot(313)
ax4.set_xlabel(r"$\phi \ [2\pi]$",fontsize=fs)
ax4.set_ylabel(r"$\delta$",fontsize=fs)
ax4.plot(r/rmax,specLoMid,'b')
ax4.plot(r/rmax,specLoHi,'r')

plt.show()

