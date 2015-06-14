#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Thu Jun 11 11:37:31 2015

@author: rafik
"""


import matplotlib as mpl
mpl.rcParams['font.family'] = 'serif'

import numpy as np
import matplotlib.pyplot as plt

from sys import argv

nplots=3

figopts = {
    'figsize': (5,3*nplots),
    'dpi': 200,
}
savefigopts = {
    'dpi': 200,
    'tightlayout':False
}
plotopts = {
    'color':'red',
    'linewidth':2
}




#def ratio_s_v(npp, nx, ny):
#	return 4.*ny/(nx/npp)/(ny+4)


if len(argv) > 2:
    print "Usage: plot_hardscaling.py data.txt"
    exit(1)
elif len(argv)==2:
    datafile = argv[1]
else:
    datafile = 'results_hard.txt'
    

data = np.loadtxt(datafile)

N = data[:,0]
t = data[:,1]


def plot_node_bounds(ax):
    for x in [2,4,8,16,32,64, 128, 256, 512]:
        ax.plot([x,x], [1e-10,1e6], 'b:')
    


fig = plt.figure(**figopts)
#fig.suptitle('Strong scaling for a 1000 x 100 grid')

### FIRST FIGURE #############################################################
ax1 = fig.add_subplot(nplots,1,1)
ax1.loglog(N, t, 'r-', **plotopts)

# plot best possible
fact=1e3
ax1.plot([1,fact], [t[0],t[0]/fact], 'k--')

# plot node boundaries
plot_node_bounds(ax1)

ax1.set_title('Strong scaling for a 1000 x 100 grid')

ax1.set_ylim(ymin=1, ymax = 1e4)
ax1.set_xlim(xmin=1, xmax=1e3)

ax1.set_xlabel('number of cpus')
ax1.set_ylabel('time [s]')


### SECOND FIGURE ############################################################
t1=t[0]

speedup = t1 / t

ax2 = fig.add_subplot(nplots,1,2)

ax2.plot(N, speedup, 'r-', **plotopts)

# plot best possible
xxyy = np.logspace(0,np.log10(fact))
ax2.plot(xxyy, xxyy, 'k--')

# plot node boundaries
plot_node_bounds(ax2)

ax2.set_xscale("log")
ax2.set_yscale("log")

ax2.set_ylim(ymin=1, ymax=1e3)
ax2.set_xlim(xmin=1, xmax=1e3)

ax2.set_title('speedup')




### 3RD FIGURE ############################################################

ax3 = fig.add_subplot(nplots,1,3)

eff = speedup / N

ax3.plot(N, eff, 'r-', **plotopts)

# plot best possible
ax3.plot([1,fact], [1,1], 'k--')

# plot node boundaries
plot_node_bounds(ax3)

ax3.set_xscale("log")
ax3.set_yscale("linear")

ax3.set_xlim(xmin=1, xmax=1e3)
ax3.set_ylim(ymin=0, ymax=1.1)

ax3.set_title('efficency')


### END ######################################################################

plt.tight_layout()

fig.savefig('hardscaling.png', **savefigopts)
plt.close(fig)




#
#
#title(r'Strong scaling for a 1000 x 100 grid', fontsize = 20)
#xlabel('Number of cores', fontsize = 20)
#ylabel('Time [s]', fontsize = 20)
#
#
#print N
#print ratio_s_v(N,1000,100)
#
#loglog(N, t, '-', color='red', linewidth=2)
#
#savefig('hardscaling.png')
#close()
#
#print "Done."
