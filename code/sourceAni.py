# animation of time-dependent source

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import DiffFunc2D as df

rmin=0.0; rmax=2*np.pi; nr=101;
Emin=-2; Emax=1; nE=100; n0=np.array(([1.0,1.0])); gamma=-0.5;

r = np.linspace(rmin,rmax,nr)
# the source distribution profile
N1 = df.initN(rmin, rmax, nr, n0[0], "twoArm")

fs=20
lw=2

fig = plt.figure(figsize=(6,3))
fig.subplots_adjust(hspace=0.45, wspace=0.4)

# First set up the figure, the axis, and the plot element we want to animate
ax = fig.add_subplot(111)
ax.set_xlim(0,1)
ax.set_ylim(0, np.max(N1))
ax.set_ylabel(r"$N$",fontsize=fs)
ax.set_xlabel(r"$\phi$",fontsize=fs)
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['left'].set_smart_bounds(True)
ax.spines['bottom'].set_smart_bounds(True)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_linewidth(lw)
ax.spines['bottom'].set_linewidth(lw)
ax.tick_params('both', length=8, width=1, which='major')
line, = ax.plot([], [], lw=lw)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    x = r/rmax
    y = np.roll(N1,-i)
    line.set_data(x, y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=2*nr, interval=100, blit=True)
anim.save("sourceAni.mp4", fps=10)

plt.show()
