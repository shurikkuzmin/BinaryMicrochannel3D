#!/usr/bin/python
import numpy
import matplotlib
matplotlib.use("GTKAgg")
import pylab
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


if __name__=="__main__":
    array=numpy.loadtxt("../phase01200.dat")
    
    ny=41
    nz=241
      
    data=array[:, 2:3]
 
    absnumpy=numpy.reshape(data,  (nz, ny))
    pylab.imshow(absnumpy)
    fig = pylab.figure()
    ax = Axes3D(fig)
    X = numpy.arange(0, ny)
    Y = numpy.arange(0, nz)
    X, Y = numpy.meshgrid(X, Y)
    ax.plot_surface(X, Y, absnumpy, rstride=1, cstride=1, cmap=cm.jet)
    #print "Velocity in the center=",absnumpy[ny/2,nx/2]

    pylab.show()
