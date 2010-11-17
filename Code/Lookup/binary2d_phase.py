#!/usr/bin/python
import numpy
import matplotlib
matplotlib.use("GTKAgg")
import pylab
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


if __name__=="__main__":
    array=numpy.loadtxt("../ZouHe/tmp/density000600.dat")
    
    ny=9
    nx=101
      
    #data=array[:, 2:3]
 
    #absnumpy=numpy.reshape(data,  (nz, ny))
    #pylab.imshow(absnumpy)
    pylab.imshow(array)
    fig = pylab.figure()
    ax = Axes3D(fig)
    X = numpy.arange(0, nx)
    Y = numpy.arange(0, ny)
    X, Y = numpy.meshgrid(X, Y)
    ax.plot_surface(X, Y,array, rstride=1, cstride=1, cmap=cm.jet)
    #print "Velocity in the center=",absnumpy[ny/2,nx/2]

    pylab.show()
