# -*- coding: utf-8 -*-
#!/usr/bin/python2

#################
# PLOT PROFILES #
#################

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
Re = 1000
y_l = 0.0
y_u = 1.0
max_step = 1

# Init stuff
X = np.zeros(1)
Y = 0.5*np.ones(1)
Xhist = []
Yhist = []
Xhist.append(X[0])
Yhist.append(Y[0])

# declare/redeclare new Disperser object from wrapped cpp, mummm cpp...
D = PyDisperser(Re, y_l, y_u)
D.gaussian_u = True
D.gaussian_v = True
D.bursting_process = False
D.bouncing_particles = False
D.direction_likelyhood = False

# update with single step until particle has reachec X=500
while (X[0] < 500.0):
    
    D.update(1, X, Y, 1.0e3, max_step)
    Xhist.append(X[0])
    Yhist.append(Y[0])

fig1, ax1 = plt.subplots(1, 1, num=1)
fig1.clf()
fig1.suptitle(r'$Re_f=$%i' % (Re))

ax1 = fig1.add_subplot(1, 1, 1)
ax1.plot(Xhist, Yhist, '-x')
# ax1.plot([0, 500], [50./Re, 50./Re], '--r')
# ax1.plot([0, 500], [100./Re, 100./Re], '--r')
ax1.set_xlim(0, 500)
ax1.set_xlabel(r'$X$')
ax1.set_ylabel(r'$Y$')

fig1.show()

fname = 'plots/ParticlePath.pdf'
plt.savefig(fname)
