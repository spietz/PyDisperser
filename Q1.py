# -*- coding: utf-8 -*-
#!/usr/bin/python2

###############
# CONVERGENCE #
###############

from mymod import PyDisperser
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.close("all")

# Customizing plots
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['font.size'] = 14

## Inputs #############################
N = 20000
Re = 1000.0
y_l = 0.0
y_u = 1.0
max_step = 100000
Um = 19.77

# New Disperser object from wrapped cpp, mummm cpp...
D = PyDisperser(Re, y_l, y_u)
D.gaussian_u = True
D.u_f = 2.5
D.gaussian_v = False
D.bursting_process = False
D.bouncing_particles = False
D.direction_likelyhood = False

# Init stuff
Y = 0.5*np.ones(N)  # use initial lump
Y = np.linspace(y_l,y_u,N)
X = np.zeros(N)

# Get results from object by taking a single step for all y
D.update(N, X, Y, 10.0, max_step)

# Calculate mean and variance
Nvec = []
Xmean = []
Xvar = []
for n in range(10, N, 100):
    vec = np.round(np.linspace(0,N-1,n).astype(int))
    Nvec.append(n)
    Xmean.append(np.mean(X[vec]))  # Xmean.append(np.mean(X[0:n]))
    Xvar.append(np.var(X[vec]))

## Plotting ###########################
fig1, axarr = plt.subplots(2, 1, sharex=True, num=1)
fig1.clf()
fig1.suptitle(r'$Re_f=$%i' % (Re))

## Xmean #######
axarr[0] = fig1.add_subplot(2, 1, 1)
axarr[0].plot(Nvec, Xmean, '-x', label=r'$\bar X_{num}$')
axarr[0].plot([1, max(Nvec)],
              [Um*10.0, Um*10.0], '--r', label=r'$U_m \Delta T$')
axarr[0].set_ylabel('$\overline X$')
axarr[0].grid(True)
axarr[0].legend(loc=0)
for label in axarr[0].get_xticklabels():
    label.set_visible(False)

## Xvar ########
axarr[1] = fig1.add_subplot(2, 1, 2)
axarr[1].plot(Nvec, Xvar, '-x')
axarr[1].set_ylabel('$\overline{(X-\overline X)^2}$')
axarr[1].set_xlabel('$N$')
axarr[1].grid(True)

fname = 'plots/convergence.pdf'
plt.savefig(fname)

fig1.show()
