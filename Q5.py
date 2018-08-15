# -*- coding: utf-8 -*-
#!/usr/bin/python2

#################
# PLOT PROFILES #
#################

from mymod import PyDisperser
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import my_fun

plt.close("all")

# Customizing plots
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['font.size'] = 14

## Inputs #############################
N = 15000
T = np.linspace(0.0, 25.0, 6)
ReVec = [1000, 5000, 10000]
y_l = 0.0
y_u = 1.0
max_step = 10000

# Init stuff
# Time between saves
dT = np.diff(T)
Xmean = np.zeros((len(dT), len(ReVec)))
Xvar = np.zeros((len(dT), len(ReVec)))
Um = np.zeros(len(ReVec))
D1 = np.zeros(len(ReVec))
TL_Xvar = np.zeros((len(dT), len(ReVec)))
TL_Xmean = np.zeros((len(dT), len(ReVec)))
Y = np.zeros((len(dT), len(ReVec), N))
X = np.zeros((len(dT), len(ReVec), N))
leg1 = []
leg2 = []

# Loop Re numbers
for r in range(0, len(ReVec)):
    Re = ReVec[r]

    # declare/redeclare new Disperser object from wrapped cpp, mummm cpp...
    D = PyDisperser(Re, y_l, y_u)
    D.gaussian_u = False
    D.u_f = False
    D.gaussian_v = False
    D.bursting_process = False
    D.bouncing_particles = True
    D.direction_likelyhood = True
    D.fall_vel = 0.15
    D.displacement = 0.0

    # set particle distribution at x=0
    Y[0, r, :] = np.linspace(y_l, y_u, N)

    # Get results for each T by updating with dT
    for t in range(0, len(dT)):

        D.update(N, X[t, r, :], Y[t, r, :], dT[t], max_step)

        if (t+1 != len(dT)):
            X[t+1, r, :] = X[t, r, :]
            Y[t+1, r, :] = Y[t, r, :]

        # Calculate mean and variance
        Xmean[t, r] = np.mean(X[t, r, :])
        Xvar[t, r] = np.var(X[t, r, :])

    # Linear regressions
    p = np.polyfit(T[1:len(T)], Xmean[:, r], 1)
    TL_Xmean[:, r] = np.polyval(p, T[1:len(T)])
    Um[r] = p[0]
    # leg1.append(r'$U_m \approx$ %0.3f for $Re_f=$%i' % (Um[r], Re))
    leg1.append(r'$Re_f=$%i' % (Re))
    
    p = np.polyfit(T[1:len(T)], Xvar[:, r], 1)
    TL_Xvar[:, r] = np.polyval(p, T[1:len(T)])
    D1[r] = 0.5*p[0]
    # leg2.append(r'$D_1  \approx$%0.3f for $Re_f=$%i' % (D1[r], Re))
    leg2.append(r'$Re_f=$%i' % (Re))



sti = 'plots/heavy'
fname = '%s/outfile' % (sti)
## Write to file #######################
f = open(fname, 'w')
f.write('Re\tUm_ana\tD1_ana\tUm_num\tD1_num\tUm_err\tD1_err\n')
for r in range(0, len(ReVec)):
    f.write('%i\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n'
            % (ReVec[r], my_fun.Um(ReVec[r]), 5.93,
               Um[r], D1[r],
               100*np.abs(my_fun.Um(ReVec[r])-Um[r])/my_fun.Um(ReVec[r]),
               100*np.abs(5.93-D1[r])/5.93))
f.close()

## Plotting ###########################
plt.ion()

## Xmean and Xvar #######
fig1, axarr1 = plt.subplots(2, 1, sharex=True, num=1)
fig1.clf()
fig1.suptitle(r'%i particles' % (N))

axarr1[0] = fig1.add_subplot(2, 1, 1)
axarr1[0].set_ylabel(r'$\overline{X}$')
axarr1[0].set_xlabel(r'$T$')
axarr1[0].grid(True)

axarr1[1] = fig1.add_subplot(2, 1, 2)
axarr1[1].set_ylabel(r'$\overline{(X-\overline{X})^2}$')
axarr1[1].set_xlabel(r'$T$')
axarr1[1].grid(True)


for r in range(0, len(ReVec)):
    axarr1[0].plot(T[1:len(T)], Xmean[:, r], 'kx')
    axarr1[0].plot(T[1:len(T)], TL_Xmean[:, r], '-', label=leg1[r])

    axarr1[1].plot(T[1:len(T)], Xvar[:, r], 'kx')
    axarr1[1].plot(T[1:len(T)], TL_Xvar[:, r], '-', label=leg2[r])

for label in axarr1[0].get_xticklabels():
    label.set_visible(False)

axarr1[0].legend(loc=0)
axarr1[1].legend(loc=0)

fname = '%s/Um_D1_N%i.pdf' % (sti,N)
plt.savefig(fname)

fig1.show()

## distributions #######
fig2, axarr2 = plt.subplots(3, 1, sharex=True, num=2)
fig2.clf()
fig2.suptitle(r'Distribution of %i particles' % (N))

for r in range(0, len(ReVec)):
    axarr2[r] = fig2.add_subplot(len(ReVec), 1, r+1, sharex=axarr2[r])
    axarr2[r].set_title(r'$Re=$%i' % (ReVec[r]))
    axarr2[r].plot(X[:, r, :].T, Y[:, r, :].T, 'x')
    axarr2[r].set_ylabel(r'$Y$')
    axarr2[r].grid(True)
    if (r != len(ReVec)-1):
        for label in axarr2[r].get_xticklabels():
            label.set_visible(False)
    axarr2[r].set_xlim(np.min(X), np.max(X))
axarr2[r].set_xlabel(r'$X$')

fname = '%s/dist_N%i.png' % (sti,N)
plt.savefig(fname)

fig2.show()


## Histograms ##
fig = plt.figure(3)
fig.clf()


# Use data for Re = 1000 and T = 10
r = 0; t = 2
n, bins, patches = plt.hist(X[t, r, :],30)
#y = mlab.normpdf(bins,Xmean[-1,r],Xvar[-1,r])
#plt.plot(bins,max(n)*y)

plt.title('Histrogram of final X position, T = %i'%T[t])
plt.ylabel('$N$')
plt.xlabel('$X$')

fname = '%s/histX.pdf' % (sti)
plt.savefig(fname)
fig.show()


fig = plt.figure(4)
fig.clf()

n, bins, patches = plt.hist(Y[t, r, :],30,orientation='horizontal')

plt.title('Histrogram of final Y position, T = %i'%T[t])
plt.xlabel('$N$')
plt.ylabel('$Y$')
fname = '%s/histY.pdf' % (sti)

plt.savefig(fname)
fig.show()


## CONCENTRATION PLOTS ##
# Using imshow
fig = plt.figure(5)
fig.clf()

H, xedges, yedges = np.histogram2d(X[t,r,:], Y[t,r,:],bins=50)
# H has to be transposed because histogram doesn't follow Cartesian convention.
# http://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram2d.html
H = H.transpose()[::-1]
np.shape(H)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
# see different interpolations for imshow here
# http://stackoverflow.com/a/14728122/1121523
plt.imshow(H, extent=extent, aspect='auto')
plt.colorbar()
plt.title('Concentration of heavy dispersant')
plt.xlabel('$X$')
plt.ylabel('$Y$')

plt.show()
fname = '%s/density.pdf' % (sti)
plt.savefig(fname)


# Using contour
fig = plt.figure(6)
fig.clf()

Xcon, Ycon = np.meshgrid(xedges[0:-1], yedges[-1:0:-1])
plt.contourf(Xcon, Ycon, H,20)
plt.colorbar()
plt.title('Concentration of heavy dispersant, T = %i'%T[t])
plt.xlabel('$X$')
plt.ylabel('$Y$')

plt.show()
fname = '%s/density_contour.pdf' % (sti)
plt.savefig(fname)
