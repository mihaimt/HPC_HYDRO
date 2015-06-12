"""
This script produces a t vs N plot for the hard scaling.
"""
#import matplotlib as mpl; mpl.rcParams['font.family'] = 'serif'
from matplotlib import *
rcParams['font.family'] = 'serif'
from matplotlib.patches import *
from numpy import *
from matplotlib.pyplot import *
from sys import *

def ratio_s_v(np, nx, ny):
	return 4.*ny/(nx/np)/(ny+4)


if len(argv) != 2:
		print "Usage: plot_hardscaling.py data.txt"
		exit(1)

datafile = argv[1]

#print "Open file: " + tipsyfile
# Load data file
data = loadtxt(datafile)

nx = data[:0]
ny = data[:1]
N = data[:,2]
t = data[:,3]




title(r'Weak scaling - ny = 100, nx = 1000*$n_{cores}$',fontsize = 20)
xlabel(r'Number of cores ($n_{cores}$)', fontsize = 20)
ylabel('Time [s]', fontsize = 20)


print N
print ratio_s_v(N,1000,100)

plot(N,t,'-',color='red',linewidth=2,label='time measurements')
vlines(8,3000,6000,color='green',linestyles='dashed', label='8 cores in a socket')
vlines(16,3000,6000,color='blue',linestyles='dashed', label='16 cores in a node (two sockets)')
ylim(3500, 6000)

#plot(N, ratio_s_v(N, 1000, 100), '-', color = 'blue', linewidth=2, label=r'$n_{comm}/n_{comp}$')




#fill_between(rhocold,ucold,color='orange')

legend(loc = 4)

savefig(datafile+'.png')
close()

print "Done."
