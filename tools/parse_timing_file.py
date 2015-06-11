#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Thu Jun 11 11:37:31 2015

@author: rafik
"""

import sys
import os
import glob

import numpy as np
import matplotlib.pyplot as plt


header = []
data = []



# plot config
width = 0.45
colors = [ # instead yellow use black... but otherwise rainbow action here
    '#ff0000', '#aaaa00', '#00ff00', '#00ffff', '#0000ff', '#ff00ff',
    '#880000', '#888800', '#008800', '#008888', '#000088', '#880088',
]


def init(fn = 'timing_0000.csv'):
    
    global data, header
    
    with open(fn) as f:
        lines = f.readlines()
        lines = [_.strip('\n') for _ in lines] # remove trailing newline
        
        # split csv data
        lines = [_.split(',') for _ in lines]
    
        # cast to float
        for i, line in enumerate(lines):
            if i==0: #header
                lines[i] = [str(_) for _ in line]
            else:
                for j, elem in enumerate(line):
                    if j==0:
                        try:
                            lines[i][j] = int(elem)
                        except ValueError:
                            lines[i][j] = str(elem)
                    elif j in [1,10,17]:
                        lines[i][j] = 0
                    else:
                        try:
                            lines[i][j] = float(elem)
                        except ValueError:
                            lines[i][j] = str(elem)
                    
            
    
    # remove header    
    header = lines[0]
    
    # remove timing point numbers
    #print header
    for i, _ in enumerate(header):
        if i not in [0,1,2,10,11,17,18]:
            header[i] = _[5:] 
    #print header
    
    data = np.array(lines[1:])





def plot_overview(steps=(1,50), imgname = 'mainloop_overview.png'):

    # ceate slice that selects the chosen steps
    sel = slice(*steps)
    
    x = data[:,0].astype(int) # stepnr
    y = [0]*10
    
    #use cols 3 - 9 (incl)
    cols = (3,10)
    for i in range(*cols):
        y[i] = data[:,i]
    bottom = y[cols[0]][sel]*0 # make an empty bar, offset for bottom start of bar
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for n, i in enumerate(range(*cols)):
        b = ax.bar(x[sel], y[i][sel], width=width, color=colors[n], bottom=bottom)
        bottom += y[i][sel]
    
    #ax.set_yscale("log")
    ax.set_xlabel('step nr')
    ax.set_ylabel('time [s]')
    
    # Put a legend to the right of the current axis
    ax.legend(header[cols[0]:cols[1]], loc='upper right')
    
    fig.savefig(imgname)
    plt.close(fig)




def plot_hydro_detail(steps=(1,9), imgname='hydro_detail.png'):

    sel = slice(*steps)
    
    x = data[:,0].astype(int) # stepnr
    
    #use cols 3 - 8 (incl)
    cols1 = (3,9)
    y1 = [0]*cols1[1]
    
    for i in range(*cols1):
        y1[i] = data[:,i]
    bottom1 = y1[cols1[0]][sel]*0 # make an empty bar, offset for bottom start of bar
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    b1 = [0] * len(range(*cols1))
    
    for n, i in enumerate(range(*cols1)):
        b1[n] = ax.bar(x[sel]-width/2, y1[i][sel], width=width, color=colors[n], bottom=bottom1)
        bottom1 += y1[i][sel]
    
    
    for cols2, k, dx in [((12,17), 6, -width),((19,24), 7, width/2.)]:
        #cols2 = (12,17)
        #cols2 = (19,24)
        y2 = [0]*(cols2[1]+1)
        
        for i in range(*cols2):
            y2[i] = data[:,i]
        bottom2 = y2[cols2[0]][sel]*0
        
        for i in range(3,k):
            bottom2 += data[:,i][sel] # make an empty bar, offset for bottom start of bar
        
        b2 = [0] * len(range(*cols2))
        
        for n, i in enumerate(range(*cols2)):
            b2[n] = ax.bar(x[sel]+dx, y2[i][sel], width=width/2, color=colors[n+6], bottom=bottom2)
            bottom2 += y2[i][sel]
    
    
    
    
    #ax.set_yscale("log")
    ax.set_xlabel('step nr')
    ax.set_ylabel('time [s]')
    
    # Put a legend to the right of the current axis
    l1 = plt.legend(b1, header[cols1[0]:cols1[1]], loc='upper right', title="mainloop")
    #ax.add_artist(l1)
    l2 = plt.legend(b2, header[cols2[0]:cols2[1]], loc='lower right', title="hydro_godunov")
    ax.add_artist(l1)
    
    fig.savefig(imgname)
    plt.close(fig)





    
def plot_log_mainloop(steps=(8,15), imgname = 'log_times_mainloop.png'):
    
    #set all values of 0 to 1e-7
    #data[data[0,:]<1e-7] = 1e-7
    
    # ceate slice that selects the chosen steps
    sel = slice(*steps)
    
    x = data[:,0].astype(int) # stepnr
    y = [0]*10
    
    #use cols 3 - 9 (incl)
    cols = (3,9)
    for i in range(*cols):
        y[i] = data[:,i]
    
    tot_width = 0.8
    ncols = cols[1]-cols[0]
    dx = 1.0*tot_width/ncols
    
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for i, c in enumerate(range(*cols)):
        b = ax.bar(x[sel]+dx*(i-ncols/2.), y[c][sel], width=dx, color=colors[i])
    
    ax.set_yscale("log")
    ax.set_ylim(ymin=1e-7)
    ax.set_xlabel('step nr')
    ax.set_ylabel('log(time) [s]')
    
    # Put a legend to the right of the current axis
    ax.legend(header[cols[0]:cols[1]], loc='upper right')
    
    fig.savefig(imgname)
    #plt.show()
    plt.close(fig)



def plot_log_hydro_godunov(steps=(8,15), imgname = 'log_times_hgodunov.png'):
    
      
    #set all values of 0 to 1e-7
    #data[data[0,:]<1e-7] = 1e-7
    
    # ceate slice that selects the chosen steps
    sel = slice(*steps)
    
    x = data[:,0].astype(int) # stepnr
    y = [0]*20
    
    #use cols 3 - 9 (incl)
    cols = (12,17)
    for i in range(*cols):
        y[i] = data[:,i]
    
    tot_width = 0.8
    ncols = cols[1]-cols[0]
    dx = 1.0*tot_width/ncols
    
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for i, c in enumerate(range(*cols)):
        b = ax.bar(x[sel]+dx*(i-ncols/2.), y[c][sel], width=dx, color=colors[i+6])
    
    ax.set_yscale("log")
    ax.set_ylim(ymin=1e-7)
    ax.set_xlabel('step nr')
    ax.set_ylabel('log(time) [s]')
    
    # Put a legend to the right of the current axis
    ax.legend(header[cols[0]:cols[1]], loc='upper right')
    
    fig.savefig(imgname)
    #plt.show()
    plt.close(fig)



def plot_time_line_mainloop_log(steps=(0,1000), imgname = 'time_line_mainloop_log.png'):
        
          
    #set all values of 0 to 1e-7
    #data[data[0,:]<1e-7] = 1e-7
    
    # ceate slice that selects the chosen steps
    sel = slice(*steps)
    
    x = data[:,0].astype(int) # stepnr
    y = [0]*20
    
    #use cols 3 - 9 (incl)
    cols = (3,9)
    for i in range(*cols):
        y[i] = data[:,i]
    
    tot_width = 0.8
    ncols = cols[1]-cols[0]
    dx = 1.0*tot_width/ncols
    
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for i, c in enumerate(range(*cols)):
        if i == 0:
            mask = y[c][sel]>=1e-6
        else:
            mask = y[c][sel]>1e-6
        xx = x[sel][mask]
        yy = y[c][sel][mask]
        b = ax.plot(xx, yy, '-', color=colors[i], )
    
    ax.set_yscale("log")
    ax.set_ylim(ymin=2e-7)
    ax.set_xlabel('step nr')
    ax.set_ylabel('log(time) [s]')
    
    # Put a legend to the right of the current axis
    ax.legend(header[cols[0]:cols[1]], loc='upper right')
    
    fig.savefig(imgname)
    #plt.show()
    plt.close(fig)


def plot_time_line_hgodunov_log(steps=(0,1000), imgname = 'time_line_hgodunov_log.png'):
        
          
    #set all values of 0 to 1e-7
    #data[data[0,:]<1e-7] = 1e-7
    
    # ceate slice that selects the chosen steps
    sel = slice(*steps)
    
    x = data[:,0].astype(int) # stepnr
    y = [0]*20
    
    #use cols 3 - 9 (incl)
    cols = (12,17)
    for i in range(*cols):
        y[i] = data[:,i]
    
    tot_width = 0.8
    ncols = cols[1]-cols[0]
    dx = 1.0*tot_width/ncols
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for i, c in enumerate(range(*cols)):
        if i == 0:
            mask = y[c][sel]>=1e-6
        else:
            mask = y[c][sel]>1e-6
        xx = x[sel][mask]
        yy = y[c][sel][mask]
        b = ax.plot(xx, yy, '-', color=colors[i+6], )
    
    ax.set_yscale("log")
    ax.set_ylim(ymin=2e-7)
    ax.set_xlabel('step nr')
    ax.set_ylabel('log(time) [s]')
    
    # Put a legend to the right of the current axis
    ax.legend(header[cols[0]:cols[1]], loc='upper right')
    
    fig.savefig(imgname)
    #plt.show()
    plt.close(fig)



def main(fn = 'timing_0000.csv', prefx=''):
    
    prefx = prefx + '_'

    init(fn)    
    
    plot_overview(steps=(1,50), imgname = prefx+'mainloop_overview.png')
    plot_hydro_detail(steps=(0,9), imgname = prefx+'hydro_detail_early.png')
    plot_hydro_detail(steps=(30,39), imgname = prefx+'hydro_detail_late.png')

    plot_log_mainloop(steps=(8,15), imgname = prefx+'log_times_mainloop.png')
    plot_log_hydro_godunov(steps=(8,15), imgname = prefx+'log_times_hgodunov.png')

    plot_time_line_mainloop_log(steps=(0,1000), imgname = prefx+'time_line_mainloop_log.png')
    plot_time_line_hgodunov_log(steps=(0,1000), imgname = prefx+'time_line_hgodunov_log.png')


def help():
    print "usage: %s timing_file.csv prefix_for_images" % sys.argv[0]


if __name__ == "__main__":
    print sys.argv
    if len(sys.argv) == 3 and os.path.isfile(sys.argv[1]):
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 2 and os.path.isfile(sys.argv[1]):
        pf = '.'.join(sys.argv[1].split('.')[:-1]) + '_'
        main(sys.argv[1], pf)
    elif len(sys.argv) == 2:
        filenames = glob.glob(sys.argv[1])
        for filename in filenames:
            print "working on:", filename
            pf = '.'.join(filename.split('.')[:-1])
            main(filename, pf)
    else:
        help()













