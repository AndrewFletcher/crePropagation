import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import DiffFunc1D as df

# params
D=1000.0
rmin=0.0; rmax=10.0; nr=101; dr=1.0; 
n0=1000; injection = "no"
nt=20000; dt=0.0001; 

# set up arrays
history = np.zeros((nr, nt))

n = df.initN(rmin, rmax, nr, n0, "triangle")
r = np.linspace(rmin, rmax, nr)
history[:,0] = n
t = np.linspace(0,nt*dt,nt)  # time

df.courantDiffusion(D, dr, dt)

i=1
while i<nt:

    k1 = dt*D*df.rhsDiff(n, nr, dr)
    nNew = n+k1/2.0
    k2 = dt*D*df.rhsDiff(nNew, nr, dr)
    n = n+k2
    
    history[:,i] = n
    
    if injection == "yes":
		n = n + df.initN(rmin, rmax, nr, n0, "triangle")
	
    i = i+1

fig = plt.figure()
fig.subplots_adjust(hspace=0.45, wspace=0.4)

ax = fig.add_subplot(3,1,1)
ax.set_ylabel("n")
ax.set_xlabel("x")
ax.plot(r, history[:,0])
ax.plot(r, history[:,(nt-1)/4])
ax.plot(r, history[:,(nt-1)/2])
ax.plot(r, history[:,3*(nt-1)/4])
ax.plot(r, history[:,(nt-1)])

ax2 = fig.add_subplot(3,2,5)
ax2.plot(t,np.sum(history, axis=0))
ax2.set_title("Total n")
ax2.set_ylabel("$\Sigma n$")
ax2.set_xlabel("t")

ax3 = fig.add_subplot(3,2,6)
ax3.plot(t,history[nr/2,:],'b-')
ax3.plot(t,history[0,:],'r-.')
ax3.plot(t,history[nr-1,:], 'r--')
ax3.set_title("n at mid-point & edges")
ax3.set_ylabel("n")
ax3.set_xlabel("t")

# Make an animation

# First set up the figure, the axis, and the plot element we want to animate
ax4 = fig.add_subplot(3,1,2)
ax4.set_xlim(rmin,rmax)
ax4.set_ylim(0, np.max(history))
ax4.set_ylabel("n")
ax4.set_xlabel("r")
line, = ax4.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    x = r
    y = history[:,i]
    line.set_data(x, y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=nt-1, interval=10, blit=True)

plt.show()
