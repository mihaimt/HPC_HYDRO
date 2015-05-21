#!/usr/bin/env python
# -*- coding: utf-8 -*-

#


import os
from sys import argv
import numpy

import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot


import sys
import numpy as np
from matplotlib import pyplot as plt
import glob
import xml.etree.ElementTree as etree
import argparse

#
# This assumes that files are numbered like this:
# outputvtk_001_000001.vts
# where the first 3 digits are the nr of the core (and thus domain)
# and the second 6 digits the step
# this is the template:
fn_tmpl  = "outputvtk_%03i_%06i.vts"
fn_tmpl2 = "outputvtk_???_%06i.vts" # this one is used for the globing only


tmpdir = '_tmp'




parser = argparse.ArgumentParser(description='Convert the output of a sim (vts files) to a movie via pngs')
parser.add_argument('-i', '--input', type=str, required=True,
                    metavar='vts_files/',
                    help='a relative path to the output folder of the simulation')
parser.add_argument('-o', '--output', type=str,required=False, default='movie.mpg',
                    metavar='movie.mpg',
                    help='the relative path to the output movie')
parser.add_argument('-n', '--nsteps', type=int, required=False, default=-1,
                    help='the number of timesteps to actually use [default = -1 = all]')
parser.add_argument('-d', '--delay', type=int, required=False, default=100,
                    help='the delay between frames in ms')
parser.add_argument('--dpi', type=int, required=False, default=100,
                    help='the resolution (dots per inch)')
parser.add_argument('--sx', type=float, required=False, default=8.0,
                    help='x size of the movie in inch')
parser.add_argument('--sy', type=float, required=False, default=6.0,
                    help='y size of the movie in inch')


args = parser.parse_args()
print os.path.abspath(args.input)

try:
    args.input = os.path.abspath(args.input)
    if not os.path.isdir(args.input):
        raise KeyError
except:
    parser.error("--input is not a valid path")




script, pathh, nsteps, delay = argv





# maps variable names in vts file to python variable / key names
valnamemap = {
    'varID':'d',
    'varIU':'u',
    'varIV':'v',
    'varIP':'p',
}

vals = valnamemap.values()

class VtsReader():
    
    def __init__(self, filename):
        self.fn = filename
        self.tree = etree.parse(filename)
        self.root = tree.getroot()
        
        # a (numpy) list of coordinates
        # coords[i] == [cx[i], cy[i]]
        # use this to assign coordinates to values
        # value at position (row) i -> is at position cx[i] / cy[i] in
        # the grid
        self.coords = None
        self.cx = None
        self.cy = None
        
        # coordinate pairs for the edges
        self.edges = None

        # the dimensions of the grids and values
        self.ex = None # the extend of the grid (#cells+1)
        self.ey = None
        
        self.nx = None # the # of cells in each dimension
        self.ny = None # also equals the number of values
        
        # and finally, the values as dict (key is the value name)
        self.values = {}


    def readDimensions(self):

        # get the number of cell edge points per dimension
        extent = map(int, self.root[0][0].attrib['Extent'].strip().split())
        self.ex = extent[1] - extent[0] + 1
        self.ey = extent[3] - extent[2] + 1
    
        # get the number of cell midpoints per dimension
        self.nx = ex-1
        self.ny = ey-1
    

    def readGrid(self):

        if not self.nX:
            self.readDimensions()

        # get the edge points
        points = self.root[0][0][0]
        coords = points[0].text.strip().split('\n')
        coords = [map(float, _.split(' ')) for _ in coords]
    
        # calculate midpoints 
        coordsMPnts = []
        cXL = []
        cYL = []
        for i in range(self.nx):
            for j in range(self.ny):
                px = (coords[i+j*ex][0] + coords[(i+1)+j*ex][0])*0.5
                py = (coords[i+j*ex][1] + coords[i+(j+1)*ex][1])*0.5
                coordsMPnts.append((px, py))
                cXL.append(px)
                cYL.append(py)
    
    
        self.edges = numpy.array(coords)
        self.coords = numpy.array(coordsMPnts)
        self.cx = numpy.array(cXL)
        self.cy = numpy.array(cYL)

    # note: for this to run, you don't need to parse the points first!!
    def readValues(self):

        celldata = root[0][0][1]
    
        for c in celldata:
            valname = c.attrib["Name"]
            valname = valnamemap.get(valname, valname)
            data = map(float, c.text.strip().split())
    
            self.values[valname] = data


'''looks at the first step and checks how many files, aka cores are around'''
def getNCores():

    lst = glob.glob(os.path.join(args.input, fn_tmpl2 % 1)) # check with step nr 1

    return len(lst)


'''looks at the first core and checks how many files, aka steps are around'''
def getNSteps():

    lst = glob.glob(os.path.join(args.input, fn_tmpl2 % 1)) # check with step nr 1

    return len(lst)



def doPlot(data, step):
    
    fig = plt.figure()
    
    for i, val in enumerate(vals):
        ax = fig.add_subplot(4,1,i)
        ax.set_title('variable %s' % val, fontsize=14)
        
        im = ax.imshow(data[val], cmap=plt.get_cmap('hot'))
        im.set_clim([0,5])
        
    fig.savefig(os.path.join(args.input, tmpdir,'img_%06i.png'%step), dpi=args.dpi)
    plt.close(fig)



def doMovie():
    pngfiles = os.path.join(args.input, tmpdir, 'img_*.png')

    exec_string = "convert %s -delay %i %s" % (pngfiles, args.delay, args.output)

    print exec_string
    os.system(exec_string)




nsteps_avail = getNSteps()
args.ncores = getNCores()

if args.nsteps <= 0 or args.nsteps > nsteps_avail:
    args.nsteps = nsteps_avail



# get the dimensions of the total grid once (expect those not to change)
# for the step nr 1
nx = 0
ny = 0
offsets = [] # at witch x=offset[i] offset start the values of core i 
subdims = [] # how many clumns have the values of core i (nx = sum_i subdims[i])
for core in range(args.ncores):
    vts = VtsReader(fn_tmpl % (core, 1))
    vts.readDimensions()
    if ny>0 and ny != vts.ny:
        raise ValueError("the files have inconsitent dimensions along the y axis!!")

    offsets.append(nx)
    subdims.append(vts.nx)
    nx = nx + vts.nx
    ny = vts.ny


# for each step (fortran counters...)
for step in range(1, args.nsteps+1):
    
    data = {}
    
    for val in vals:
        data[val] = np.zeros((0, ny))
    
    for core in range(args.ncores):
        vts = VtsReader(fn_tmpl % (core, step))
        vts.readValues()
        
        for val in vals:
            tmp = vts.values[val].reshape((subdims[core], ny))
            data[val] = np.hstack((data[val], tmp))
            
    doPlot(data, step, val)


doMovie()


