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

if len(argv) != 2:
		print "Usage: plot_hardscaling.py data.txt"
		exit(1)

datafile = argv[1]

#print "Open file: " + tipsyfile
# Load data file
data = loadtxt(datafile)

N = data[:,0]
t = data[:,1]

title(r'Hard scaling for a 1000 x 100 grid')
xlabel('Number of cores')
ylabel('Time [s]')

scl = 100.
loglog(N,t,'-',color='red',linewidth=1,label='zbox')
loglog([1,scl], [t[0], t[0]/scl],':',color='blue',linewidth=1)

#fill_between(rhocold,ucold,color='orange')

legend()

savefig(datafile+'.png')
close()

print "Done."
