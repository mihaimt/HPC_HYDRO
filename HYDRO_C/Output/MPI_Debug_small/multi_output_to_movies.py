import os
from sys import argv
import numpy
from matplotlib import pyplot

def add_zeros(step):
	
	if int(step) < 10:
		outt = '0000'+str(step)
	elif int(step) < 100:
		outt = '000' +str(step)
	elif int(step) < 1000:
		outt = '00' + str(step)
	elif int(step) < 10000:
		outt = '0'  + str(step)
	else:
		outt = str(step)
	
	return outt

print "="*40
print "proper call:"
print "python multi_output_to_movies.py <path to output directory> <number of steps> <delay between movie frames>"
print "-"*40
print "example:"
print "python multi_output_to_movies.py ~/git/HPC_Hydro/HYDRO_C/Output/Input_sedov_100x10 10 100"
print "="*40


script, pathh, nsteps, delay = argv




import xml.etree.ElementTree as etree

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

    values['coordsEdges'] = numpy.array(coords)
    values['coordsMidPoints'] = numpy.array(coordsMPnts)
    values['coordsMidPointsX'] = numpy.array(cXL)
    values['coordsMidPointsY'] = numpy.array(cYL)

    for c in celldata:
        valname = c.attrib["Name"]
        valname = valnamemap.get(valname, valname)
        data = map(float, c.text.strip().split())

        values[valname] = data
    return values, nx, ny, filename





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
#	print large_d	

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
	 
	#print large_d



