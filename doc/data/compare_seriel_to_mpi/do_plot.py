#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 19:34:06 2015

@author: rafik
"""

import numpy as np
import matplotlib.pyplot as plt


dfile = 'data.txt'


with open(dfile) as f:
    lns = f.readlines()
    
times = []

for ln in lns:
    try:
        t = float(ln.split('(')[1].split(')')[0])
    except IndexError:
        continue
    times.append(t)

# orginal code
t1 = np.array(times[10:15])
# orginal code, without fprint
t2 = np.array(times[15:20])
# new code USE_MPI OFF
t3 = np.array(times[0:5])
# new code USE_MPI ON
t4 = np.array(times[5:10])


tt = [t1,t2,t3,t4]

tm = [0,]*4
te = [0,]*4
for i, t in enumerate(tt):
    print i, t
    tm[i] = np.mean(t)
    te[i] = np.std(t)

tm = np.array(tm)
te = np.array(te)

rel = tm / tm[0]
relerr = te / tm[0]

fig = plt.figure(figsize=(4.7,3), dpi=200)
ax = fig.add_subplot(111)

x = np.arange(len(tt)) + 1

width = 0.8

ax.bar(x-width/2., rel, width=width)

ax.set_xticks(x)
ax.set_xticklabels( ('mono', 'mono\nno fprintf', 'parallel\nUSE_MPI 0', 'parallel\nUSE_MPI 1',) )
ax.set_ylim(ymax=1.05, ymin=0.8)
ax.set_ylabel('relative runtime [1]')
ax.set_title('Comparison of codes')

plt.tight_layout()

fig.savefig('plot.svg', dpi= 200)
fig.savefig('plot.png', dpi= 200)
