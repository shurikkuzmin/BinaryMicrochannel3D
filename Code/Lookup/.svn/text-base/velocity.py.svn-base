#!/usr/bin/python
import numpy
import matplotlib
matplotlib.use("GTKAgg")
import pylab
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


if __name__=="__main__":
    array=numpy.loadtxt("../phase00600.dat")
    
    nx=3
    ny=81
    nz=241
    
    slice=nz/2
    
    data=array[slice*nx*ny:(slice+1)*nx*ny, 3:6]
    absmod=[]
    for counter in range(0, nx*ny):
        absmod.append(numpy.sqrt(data[counter,0]**2+data[counter, 1]**2+data[counter, 2]**2))
    absnumpy=numpy.array(absmod)
    absnumpy=numpy.reshape(absnumpy, (ny, nx))
    pylab.imshow(absnumpy)
    fig = pylab.figure()
    ax = Axes3D(fig)
    X = numpy.arange(0, nx)
    Y = numpy.arange(0, ny)
    X, Y = numpy.meshgrid(X, Y)
    ax.plot_surface(X, Y, absnumpy, rstride=1, cstride=1, cmap=cm.jet)
    print "Velocity in the center=",absnumpy[ny/2,nx/2]

    pylab.show()
