# -*- coding: utf-8 -*-
#!/usr/bin/python2


# Animation related
# http://matplotlib.org/api/animation_api.html
# http://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/
# http://stackoverflow.com/questions/17212722/matplotlib-imshow-how-to-animate
# http://matplotlib.org/examples/animation/

#################
# PLOT PROFILES #
#################

from mymod import PyDisperser
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import colors

plt.close("all")


## Inputs #############################
N = 1000
T = np.linspace(0.0, 10.0, 300)
T = np.linspace(0.0, 5.0, 150)
# Remember that beta<1. beta = W/kappa, kappa = 0.41 (eq. 4.102)
WVec = np.linspace(0, 0.25, 15)
WVec = [0.15]
y_l = 0.0
y_u = 1.0
max_step = 1500
Re = 1000

# Init stuff
# Time between saves
dT = np.diff(T)
Xmean = np.zeros((len(dT), len(WVec)))
Xvar = np.zeros((len(dT), len(WVec)))
Um = np.zeros(len(WVec))
D1 = np.zeros(len(WVec))
Y = np.zeros(N)
X = np.zeros(N)

# In order to use the same color scale for all plots, we need this
# http://matplotlib.org/examples/pylab_examples/multi_image.html
class ImageFollower:
    'update image in response to changes in clim or cmap on another image'
    def __init__(self, follower):
        self.follower = follower
    def __call__(self, leader):
        self.follower.set_cmap(leader.get_cmap())
        self.follower.set_clim(leader.get_clim())

# Loop Re numbers
for r in range(0, len(WVec)):
    W = WVec[r]

    # declare/redeclare new Disperser object from wrapped cpp, mummm cpp...
    D = PyDisperser(Re, y_l, y_u)
    D.gaussian_u = False
    D.u_f = False
    D.gaussian_v = False
    D.bursting_process = False
    D.bouncing_particles = True
    D.direction_likelyhood = True
    D.fall_vel = W
    D.displacement = 0.0

    # set particle distribution at x=0
    Y = np.linspace(y_l, y_u, N)

    # ims is a list of lists, each row is a list of artists to draw in the
    # current frame; here we are just animating one artist, the image, in
    # each frame
    ims = []
    extent = [0, 340, 0, 1]
    vmin = 1e40
    vmax = -1e40
    fig = plt.figure()
    # Get results for each T by updating with dT
    for t in range(0, len(dT)):

        D.update(N, X, Y, dT[t], max_step)

        # Calculate mean and variance
        Xmean[t, r] = np.mean(X)
        Xvar[t, r] = np.var(X)

        ## plot
        # histogram needs the same range as imshow should plot.
        y_dist = np.diff(extent[2:])
        x_dist = np.diff(extent[:2])
        H, xedges, yedges = np.histogram2d(X, Y, bins=[x_dist[0]*20, y_dist[0]*20],
                                           range=[extent[0:2], extent[2:4]])
        H = H.transpose()[::-1]
        dd = np.ravel(H)
        # Manually find the min and max of all colors for
        # use in setting the color scale.
        vmin = min(vmin, np.amin(dd))
        vmax = max(vmax, np.amax(dd))
        im = plt.imshow(H, extent=extent, aspect='auto')
        ims.append([im])

    vmax = 10
    ## Create animation
    # Set the color scale
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    for i, im in enumerate(ims):
        im[0].set_norm(norm)
        if i > 0:
            ims[0][0].callbacksSM.connect('changed', ImageFollower(im[0]))

    plt.colorbar(ims[0][0])
    ani = animation.ArtistAnimation(fig, ims, blit=False)
    # ani.save('dispersion_W=%0.3f.mp4'%(W) , fps=30,
    #           extra_args=['-vcodec', 'libx264'])
    plt.show()
