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

nplots=2
fact=1e2


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
    print "Usage: %s data.txt" % argv[0]
    exit(1)
elif len(argv)==2:
    datafile = argv[1]
else:
    datafile = 'results_weak.txt'
    

data = np.loadtxt(datafile)

nx = data[:0]
ny = data[:1]
N = data[:,2]
t = data[:,3]


def plot_node_bounds(ax):
    for x in [2,4,8,16,32,64, 128, 256, 512]:
        ax.plot([x,x], [1e-10,1e6], 'b:')
    


fig = plt.figure(**figopts)


### FIRST FIGURE #############################################################
ax1 = fig.add_subplot(nplots,1,1)

ax1.semilogx(N, t, 'r-', **plotopts)

# plot best possible
ax1.plot([1,fact],[t[0],t[0]], 'k--')

# plot node boundaries
plot_node_bounds(ax1)

ax1.set_title('Weak scaling ($n_x=1000 \cdot n_{cores}$; $n_y=100$)')

ax1.set_xlim(xmin=1, xmax=fact)
ax1.set_ylim(ymin=3000, ymax = 6000)

ax1.set_xlabel('number of cpus')
ax1.set_ylabel('time [s]')




### SECOND FIGURE ############################################################
t1=t[0]

eff = t1 / t

ax2 = fig.add_subplot(nplots,1,2)

ax2.plot(N, eff, 'r-', **plotopts)

# plot best possible
ax2.plot([1,fact], [1,1], 'k--')

# plot node boundaries
plot_node_bounds(ax2)

ax2.set_xscale("log")
ax2.set_yscale("linear")

ax2.set_xlim(xmin=1, xmax=fact)
ax2.set_ylim(ymin=0, ymax=1.1)

ax2.set_title('efficency')




#### 3RD FIGURE ############################################################
#
#ax3 = fig.add_subplot(nplots,1,3)
#
#eff = speedup / N
#
#ax3.plot(N, eff, 'r-', **plotopts)
#
## plot best possible
#ax3.plot([1,fact], [1,1], 'k--')
#
## plot node boundaries
#plot_node_bounds(ax3)
#
#ax3.set_xscale("log")
#ax3.set_yscale("linear")
#
#ax3.set_xlim(xmin=1, xmax=1e3)
#ax3.set_ylim(ymin=0, ymax=1.1)
#
#ax3.set_title('efficency')





### END ######################################################################

plt.tight_layout()

fig.savefig('weakscaling.svg', **savefigopts)
fig.savefig('weakscaling.png', **savefigopts)
plt.close(fig)

