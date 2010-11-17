#/usr/bin/python

import numpy
from enthought.mayavi.mlab import *

def test_quiver3d():
    #x, y, z = numpy.mgrid[-1:1, -1:1, -1:1]
    x=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    y=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    z=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    u=numpy.array([0,1,-1,0, 0,0, 0,1,-1, 1,-1,0, 0, 0, 0,1,-1, 1,-1])
    v =numpy.array([0,0, 0,1,-1,0, 0,1, 1,-1,-1,1,-1, 1,-1,0, 0, 0, 0])
    w=numpy.array([0,0, 0,0, 0,1,-1,0, 0, 0, 0,1, 1,-1,-1,1, 1,-1,-1])

    

    obj = quiver3d(x, y, z, u, v, w, line_width=3, scale_factor=1)
    for counter in range(0, len(u)):
        text3d(u[counter], v[counter], w[counter], str(counter), scale=0.1)
    axes(extent=[-1, 1, -1, 1, -1,1], xlabel="X", ylabel="Y", zlabel="Z")
    return obj

test_quiver3d()
show()
#import matplotlib as mpl
#from mpl_toolkits.mplot3d import Axes3D
#import numpy as np
#import matplotlib.pyplot as plt
#
#mpl.rcParams['legend.fontsize'] = 10
#
#fig = plt.figure()
#ax = Axes3D(fig)
#theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
##z = np.linspace(-2, 2, 100)
##r = z**2 + 1
##x = r * np.sin(theta)
##y = r * np.cos(theta)
#x=[0, 1, 0, 1]
#y=[0, 1, 0, -1]
#z=[0, 1, 0, 0]
#ax.plot(x, y, z, label="parametric curve")
#
#ax.set_xlabel("X Label")
#ax.set_ylabel("Y Label")
#ax.set_zlabel("Z Label")
#
#plt.show()
