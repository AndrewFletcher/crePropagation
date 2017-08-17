# code to make some figures for talks

import numpy as np
import matplotlib.pyplot as plt
import DiffFunc2D as df

rmin=0.0; rmax=2*np.pi; nr=101;
Emin=-2; Emax=1; nE=100; n0=np.array(([1.0,1.0])); gamma=-0.5;

x = np.linspace(rmin,rmax,nr)
E = np.logspace(Emin,Emax,nE)
# the source distribution profile
N1 = df.initN(rmin, rmax, nr, n0[0], "twoArm")
# the energy spectrum
N2=df.initE(E, n0, gamma)[0,:]

fs=20
lw=2

fig=plt.figure(figsize=(6,8))
fig.subplots_adjust(hspace=0.45, wspace=0.4)

ax1=fig.add_subplot(211)
ax1.plot(x/(rmax),N1,linewidth=lw)
ax1.set_xlabel(r'$\phi$',fontsize=fs)
ax1.set_ylabel(r'$N$',fontsize=fs)
ax1.spines['right'].set_color('none')
ax1.spines['top'].set_color('none')
ax1.spines['left'].set_smart_bounds(True)
ax1.spines['bottom'].set_smart_bounds(True)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
ax1.spines['left'].set_linewidth(lw)
ax1.spines['bottom'].set_linewidth(lw)
ax1.tick_params('both', length=8, width=1, which='major')

ax2=fig.add_subplot(212)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.plot(E,N2,linewidth=lw)
ax2.set_xlabel(r'$E$',fontsize=fs)
ax2.set_ylabel(r'$N$',fontsize=fs)
ax2.spines['right'].set_color('none')
ax2.spines['top'].set_color('none')
ax2.spines['left'].set_smart_bounds(True)
ax2.spines['bottom'].set_smart_bounds(True)
ax2.xaxis.set_ticks_position('bottom')
ax2.yaxis.set_ticks_position('left')
ax2.spines['left'].set_linewidth(lw)
ax2.spines['bottom'].set_linewidth(lw)
ax2.tick_params('both', length=10, width=1, which='major')
ax2.tick_params('both', length=5, width=1, which='minor')

plt.savefig('source.jpg')

plt.show()



