#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Achtung, die koordinaten sind cell point coordinates

die anderen values sind values im innern einer zelle!

"""

import sys
import xml.etree.ElementTree as etree  
import numpy as np

# maps variable names in vts file to python variable / key names
valnamemap = {
    'varID':'d',
    'varIU':'u',
    'varIV':'v',
    'varIP':'p',
}


def read_vtk(filename):

    values = {}
    
    tree = etree.parse(filename)
    root = tree.getroot()

    # get the number of cell edge points per dimension
    extent = map(int, root[0][0].attrib['Extent'].strip().split())
    ex = extent[1] - extent[0] + 1
    ey = extent[3] - extent[2] + 1
    
    # get the number of cell midpoints per dimension
    nx = ex-1
    ny = ey-1
    
    # get the edge points
    points = root[0][0][0]
    coords = points[0].text.strip().split('\n')
    coords = [map(float, _.split(' ')) for _ in coords]
    
    # calculate midpoints 
    coordsMPnts = []
    cXL = []
    cYL = []
    for i in range(nx):
        for j in range(ny):
            px = (coords[i+j*ex][0] + coords[(i+1)+j*ex][0])*0.5
            py = (coords[i+j*ex][1] + coords[i+(j+1)*ex][1])*0.5
            coordsMPnts.append((px, py))
            cXL.append(px)
            cYL.append(py)
    
    celldata = root[0][0][1]
    
    values['coordsEdges'] = np.array(coords)
    values['coordsMidPoints'] = np.array(coordsMPnts)
    values['coordsMidPointsX'] = np.array(cXL)
    values['coordsMidPointsY'] = np.array(cYL)

    for c in celldata:
        valname = c.attrib["Name"]
        valname = valnamemap.get(valname, valname)
        data = map(float, c.text.strip().split())
        
        values[valname] = data
    return values



def read_many_vtk(filenames):
    data = []
    for fn in filenames:
        d = read_vtk(fn)
        data.append(d)
    return data



if __name__ == "__main__":
    args = sys.argv
    if len(args)<2:
        print "\nusage: %s filename1.vtk [filename2.vtk, [...]]\n" % args[0]
    else:
        files = read_many_vtk(sys.argv[1:])
        for f in files:
            
            print "modify stdoutput here"
            
            for i, coords in enumerate(f['coordsMidPoints']):
                print i, coords, f['d'][i], f['p'][i], f['u'][i], f['v'][i]
