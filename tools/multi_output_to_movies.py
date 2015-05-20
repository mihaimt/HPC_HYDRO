#!/usr/bin/env python
# -*- coding: utf-8 -*-

#


import os
from sys import argv
import numpy
from matplotlib import pyplot


import sys
import numpy as np
from matplotlib import pyplot as plt
import glob
import xml.etree.ElementTree as etree




print "="*40
print "proper call:"
print "python multi_output_to_movies.py <path to output directory> <number of steps> <delay between movie frames>"
print "-"*40
print "example:"
print "python multi_output_to_movies.py ~/git/HPC_Hydro/HYDRO_C/Output/Input_sedov_100x10 10 100"
print "="*40



script, pathh, nsteps, delay = argv





# maps variable names in vts file to python variable / key names
valnamemap = {
    'varID':'d',
    'varIU':'u',
    'varIV':'v',
    'varIP':'p',
}

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

'''
looks at the first step and checks how many files, aka cores are around
'''
def getNSteps(dirr):

    lst = glob.glob(os.path.join(dirr, fn_tmpl2 % 1)) # check with step nr 1

    return len(lst)



#
# This assumes that files are numbered like this:
# outputvtk_001_000001.vts
# where the first 3 digits are the nr of the core (and thus domain)
# and the second 6 digits the step
# this is the template:
fn_tmpl  = "outputvtk_%03i_%06i.vts"
fn_tmpl2 = "outputvtk_???_%06i.vts" # this one is used for the globing only



nsteps = getNSteps(d)








dir_list = os.listdir(str(pathh))

step_list = range(int(nsteps)+1)


if pathh[-1] == '/':
    ppath = pathh
else:
    ppath = pathh+'/'

print "-"*40
print "Searching directory for numbered outputs: ", str(pathh)
print "-"*40

print dir_list

step_f = [[]]*(int(nsteps)+1)

i = 1

for step in step_list:
    step_string = add_zeros(step)
    
    for element in dir_list:
        if (element.find(step_string) > 0) and (element.find('.png') < 0) and (element.find('movie') < 0):
            step_f[i-1] = step_f[i-1] + [element]
    i = i +1

print "-"*40
print "Grouping outputs by step number"
print "-"*40

for element in step_f[1:]:
    element.sort()
    print element


for element in step_f[1:]:
    element.sort()
    large_d, large_p, large_u, large_v, vtk_name = [], [], [], [], []
 
    for vtk in element:

            f, nx, ny, filename = read_vtk(ppath+vtk)
            d=numpy.zeros((nx, ny))
            p=numpy.zeros((nx, ny))
            u=numpy.zeros((nx, ny))
            v=numpy.zeros((nx, ny))
            x_r = numpy.zeros((nx))
            y_r = numpy.zeros((ny))
            for i, coords in enumerate(f['coordsMidPoints']):

                    d[i%nx,i/nx] = f['d'][i]
                    p[i%nx,i/nx] = f['p'][i]
                    u[i%nx,i/nx] = f['u'][i]
                    v[i%nx,i/nx] = f['v'][i]

                    if i<nx:
                            x_r[i] = coords[0]

                    if i%nx == 0:
                            y_r[i/nx] = coords[1]
    
        large_d = large_d + [d.T] 
        large_p = large_p + [p.T]
        large_u = large_u + [u.T]
        large_v = large_v + [v.T]
        vtk_name = vtk_name + [vtk]    
        outname = ppath+vtk[0:-6]

    da = large_d[0]
#    print large_d    

    for element in large_d[1:]:
    
        dat = numpy.hstack((da, element))
        da = dat


    da = numpy.matrix(da)
    pyplot.imshow(da)
    pyplot.hot()    
    pyplot.clim(0,5)
    pyplot.title("d")
    pyplot.colorbar(orientation = "horizontal")
    pyplot.savefig(outname +"_d"+".png")
    pyplot.cla()
    pyplot.clf()
     


    pa = large_p[0]
        for element in large_p[1:]:

                pat = numpy.hstack((pa, element))
                pa = pat


        pa = numpy.matrix(pa)
        pyplot.imshow(pa)
        pyplot.hot()
        pyplot.clim(0,400)
    pyplot.title("p")    
    pyplot.colorbar(orientation = "horizontal")
        pyplot.savefig(outname +"_p"+".png")
        pyplot.cla()
        pyplot.clf()


    ua = large_u[0]
        for element in large_u[1:]:

                uat = numpy.hstack((ua, element))
                ua = uat


        ua = numpy.matrix(ua)
        pyplot.imshow(ua)
        pyplot.hot()
        pyplot.clim(0,40)
    pyplot.title("u")
    pyplot.colorbar(orientation = "horizontal")
        pyplot.savefig(outname +"_u"+".png")
        pyplot.cla()
        pyplot.clf()



        va = large_v[0]
        for element in large_v[1:]:

                vat = numpy.hstack((va, element))
                va = vat


        va = numpy.matrix(va)
        pyplot.imshow(va)
        pyplot.hot()
        pyplot.clim(-10,40)
    pyplot.colorbar(orientation = "horizontal")
    pyplot.title("v")
        pyplot.savefig(outname +"_v"+".png")
        pyplot.cla()
        pyplot.clf()




    #print large_d


print "-"*40
print "Created plots"
print "-"*40


print "-"*40
print "Creating movies"
print "-"*40


for q in ["d", "p", "u", "v"]:
    exec_string = "convert " + ppath+"*"+q+".png -delay "+str(int(delay)) +" "+ ppath+"movie_"+q+".mpg"
    print exec_string
    os.system(exec_string)
    del exec_string



print "-"*40
print "Done"
print "-"*40

print "="*40
print "The files have been placed in the directory: "
print  ppath
print "="*40

