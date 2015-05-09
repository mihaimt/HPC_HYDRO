import sys
import xml.etree.ElementTree as etree  
import numpy as np
from matplotlib import pyplot

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
    print filename
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
    return values, nx, ny, filename



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
        f, nx, ny, filename = read_vtk(sys.argv[1])
        d=np.zeros((nx, ny))
       	p=np.zeros((nx, ny))
	u=np.zeros((nx, ny))
	v=np.zeros((nx, ny))      
	x_r = np.zeros((nx))
	y_r = np.zeros((ny))
        for i, coords in enumerate(f['coordsMidPoints']):
        #        print i, coords, f['d'][i], f['p'][i], f['u'][i], f['v'][i]

		d[i%nx,i/nx] = f['d'][i]
		p[i%nx,i/nx] = f['p'][i]
		u[i%nx,i/nx] = f['u'][i]
		v[i%nx,i/nx] = f['v'][i]

		if i<nx:
			x_r[i] = coords[0]
		
		if i%nx == 0:
			y_r[i/nx] = coords[1]

#		Density
	pyplot.imshow(d.T)
	pyplot.hot()
	pyplot.title("d")
	pyplot.savefig(str(filename)+"_"+"d"+".png") 
	pyplot.cla()
	pyplot.clf()
		
#		Pressure
	pyplot.imshow(p.T)
        pyplot.hot()
	pyplot.title("p")
        pyplot.savefig(str(filename)+"_"+"p"+".png") 
        pyplot.cla()
        pyplot.clf()
	
	pyplot.imshow(u.T)
        pyplot.hot()
	pyplot.title("u")
        pyplot.savefig(str(filename)+"_"+"u"+".png") 
        pyplot.cla()
        pyplot.clf()
		
	pyplot.imshow(v.T)
        pyplot.hot()
	pyplot.title("v")
        pyplot.savefig(str(filename)+"_"+"v"+".png") 
        pyplot.cla()
        pyplot.clf()




