# energy losses test cases

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import DiffFunc1D as df

Emin=-2; Emax=2; nE = 100; n0 = 1; gamma = -2.5; 
dt = 0.00001; nt = 100

Phi = 10.0
injection='no'

E = np.logspace(Emin,Emax,nE)  # includes the ghost zone, E_min
N = df.initE(E[1:], n0, gamma)  # excludes the ghost zone N(E_min)
deltaE = E[1:]-E[:-1]  # the delta Es are [E2-E1,E3-E2,...]
df.courantEnergy(Phi, dt, deltaE)

# set up arrays
t = np.linspace(0,nt*dt,nt)  # time
history = np.zeros((nE-1,nt))  # nE-1 as nE[0] is ghost
history[:,0] = N  # initial distribution

i=1
while i<nt:
    k1 = dt*df.rhsE(Phi,N,E,deltaE)
    nNew = N+k1/2.0
    k2 = dt*df.rhsE(Phi,nNew,E,deltaE)
    N = N+k2
    N = df.floorN(N)
    history[:,i] = N
    
    if injection == "yes":
		N = N + df.initE(E[1:],n0,gamma)
	
    i = i+1

# calculate the spectral index
specLoMid = df.specInd(history[0,:],history[nE/2,:],E[1],E[1+nE/2])
specLoHi = df.specInd(history[0,:],history[nE-2,:],E[1],E[nE-1])
specMidHi = df.specInd(history[nE/2,:],history[nE-2,:],E[1+nE/2],E[nE-1])
specHi = df.specInd(history[nE-3,:],history[nE-2,:],E[nE-2],E[nE-1])

fig = plt.figure()
fig.subplots_adjust(hspace=0.45, wspace=0.4)

ax = fig.add_subplot(3,1,1)
ax.set_ylabel("n")
ax.set_xlabel("E")
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(E[1:], history[:,0])
ax.plot(E[1:], history[:,(nt-1)/4])
ax.plot(E[1:], history[:,(nt-1)/2])
ax.plot(E[1:], history[:,3*(nt-1)/4])
ax.plot(E[1:], history[:,(nt-1)])

ax2 = fig.add_subplot(3,2,5)
ax2.set_yscale('log')
ax2.plot(t,np.sum(history, axis=0))
ax2.set_title("Total n")
ax2.set_ylabel("$\Sigma n$")
ax2.set_xlabel("t")

ax3 = fig.add_subplot(3,2,6)
ax3.plot(t,history[0,:],'r-')
ax3.plot(t,history[nE/2,:], 'r-.')
ax3.plot(t,history[nE-2,:], 'r--')
ax3.set_title("n at max-E, mid-E min-E")
ax3.set_ylabel("n")
ax3.set_xlabel("t")

ax4=fig.add_subplot(3,1,2)
ax4.plot(t,specLoHi,'b')
ax4.plot(t,specLoMid,'r')
ax4.plot(t,specHi,'k')
ax4.set_ylabel(r"$\delta$")
ax4.set_xlabel(r"$t$")
ax4.set_ylim(0,5)

# Make an animation

# First set up the figure, the axis, and the plot element we want to animate
#~ ax4 = fig.add_subplot(3,1,2)
#~ ax4.set_xlim(Emin,Emax)
#~ ax4.set_ylim(0, np.max(history))
#~ ax4.set_xscale('log')
#~ ax4.set_yscale('log')
#~ ax4.set_ylabel("n")
#~ ax4.set_xlabel("E")
#~ line, = ax4.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    x = E[1:]
    y = history[:,i]
    line.set_data(x, y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
#~ anim = animation.FuncAnimation(fig, animate, init_func=init,
                               #~ frames=nt-1, interval=10, blit=True)

plt.show()

