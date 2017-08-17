# code for simple modelling of cr electron propagation
import numpy as np
import matplotlib.pyplot as plt
import DiffFunc2D as df

def newE(E,t,beta):
    # forward energy loss
    #   energy at time t, given E(t=0)
    return E/(1+beta*E*t)
    
def origE(E,t,beta):
    # reverse energy loss
    #   energy of E at time 0
    return E/(1-beta*E*t)
    
def N(E,gamma0,K):
    # number of cre at energy E
    return  K*E**(-gamma0)

def spxE(N1,N2,E1,E2):
    # the energy spectral index
    return -np.log(N1/N2)/np.log(E1/E2)
    
# params
K=1.0  # N(E)=KE**gamma0
gamma0=2.0  # spectral index at source
L=1.0  # distance, source to target
D=4.0  # diffusion coefficient
beta=1.0  # energy loss rate
E1=1.0  # energy at source
E2=4.0  # energy at source
t=L**2/(2*D)  # time for diffusion from source to target

# source, Q
N1=N(E1,gamma0,K)  # N at E1 at source
N2=N(E2,gamma0,K)  # N at E2 at source

#-------------------------------
# using the approximate formulae
# target, A
E1A=newE(E1,t,beta)  # cre with E1 at source have E1A at target
E2A=newE(E2,t,beta)
N1A=N1  # N(E1A, A) = N(E1,Q)
N2A=N2
gammaA=spxE(N1A,N2A,E1A,E2A)
E1Q=origE(E1,t,beta)  # cre with E1 at target had E1Q at source
E2Q=origE(E2,t,beta)
N1Q=N(E1Q,gamma0,K)  # N(E1,A) = N(E1Q,Q)
N2Q=N(E2Q,gamma0,K)
gammaB=spxE(N1Q,N2Q,E1,E2)

print 'gammaA ', gammaA
print 'gammaA/gamma0 ', gammaA/gamma0
print 'gammaA/gamma0 eqn. ', np.log(E1/E2)/np.log(2*D*E1/E2)
print 'gammaB ', gammaB
print 'gammaB/gamma0 ', gammaB/gamma0

#--------------------
# numerical solutions
# params
rmin=0.0; rmax1=L; nr=101; dr=L/nr;
n0 = K; 
Emin=-1; Emax1=8; nE = 100; gamma = -gamma0; 
Phi = 1.0
nt=100; dt=0.00001; 
#escT = 0.003  # escape time in time units (same as dt)
source = "yes"
sourceTime = 0.0  # how many times source will pass through domain
record = "no"

# set up arrays
history = np.zeros((nr, nE-1, nt))  # nE-1 as nE[0] is ghost
r = np.linspace(rmin, rmax1, nr)
E = np.logspace(Emin,Emax1,nE)  # includes the ghost zone, E_min
deltaE = E[1:]-E[:-1]  # the delta Es are [E2-E1,E3-E2,...]
t = np.linspace(0,nt*dt,nt)  # time
nTotal = np.zeros((nt))  # store total number of CRE
nEmax1 = np.zeros((nt)) # store N(:,Emax1)
nEmid = np.zeros((nt)) #       N(:,Emid)
nEmin = np.zeros((nt)) #       N(:,Emin)

# test stability
df.courantDiffusion(D, dr, dt)
df.courantEnergy(Phi, dt, deltaE)

# source distribution
NN0 = df.initN(rmin, rmax1, nr, n0, "delta")
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
    
 #   n = n*(1.0-(dt/escT))  # escape losses
 
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

# plotting params
fs=20
lw=2

fig=plt.figure()

ax1=fig.add_subplot(211)
ax1.plot([E1,E2,E1Q,E2Q],[N1,N2,N1Q,N2Q],'-or',label='source')
ax1.plot([E1A,E2A],[N1A,N2A],'-ob',label='target: same N')
ax1.plot([E1Q,E2Q],[N1Q,N2Q],'og')
ax1.plot([E1,E2],[N1Q,N2Q],'-og',label='target: same E')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim([E1A/2,2*E2Q])
ax1.set_ylim([N2Q/2,2*N1])
ax1.set_xlabel('Energy')
ax1.set_ylabel('Number')
ax1.legend()

ax2 = fig.add_subplot(212)
ax2.set_title("Numerics")
ax2.set_xlabel(r"$E$",fontsize=fs)
ax2.set_ylabel(r"$N$",fontsize=fs)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.plot(E[1:],n[nr-1,:],'b',linewidth=lw)

plt.show()
