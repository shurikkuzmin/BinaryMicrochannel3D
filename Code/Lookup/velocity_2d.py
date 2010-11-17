#!/usr/bin/python
import numpy
import matplotlib
matplotlib.use("GTKAgg")
import pylab
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


if __name__=="__main__":
    array=numpy.loadtxt("../phase00200.dat")
    
    ny=41
    nz=241
      
    data=array[:, 2:5]
 
    absmod=[]
    for counter in range(0, ny*nz):
        absmod.append(numpy.sqrt(data[counter,0]**2+data[counter, 1]**2+data[counter, 2]**2))
    absnumpy=numpy.array(absmod)
    absnumpy=numpy.reshape(absnumpy,  (nz, ny))
    pylab.imshow(absnumpy)
    fig = pylab.figure()
    ax = Axes3D(fig)
    X = numpy.arange(0, ny)
    Y = numpy.arange(0, nz)
    X, Y = numpy.meshgrid(X, Y)
    ax.plot_surface(Y, X, absnumpy, rstride=1, cstride=1, cmap=cm.jet)
    #print "Velocity in the center=",absnumpy[ny/2,nx/2]

    pylab.show()
