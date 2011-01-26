#!/usr/bin/python

def read(name):
    import numpy
    import pylab
    
    arr=numpy.load(name)
    
    dims=arr['phi'].shape
    velx=arr['v'][0,:,:]
    vely=arr['v'][1,:,:]
    
    print dims
    
    x,y=numpy.mgrid[0:dims[0],0:dims[1]]
    x_short=x[::30,::150]
    y_short=y[::30,::150]
    velx_short=velx[::30,::150]
    vely_short=vely[::30,::150]
  
    pylab.figure()
    pylab.imshow(arr['phi'])
    
    pylab.figure(figsize=[31,2.5])
    pylab.quiver(y_short,x_short,velx_short,vely_short,scale=0.1)
    
    pylab.show()

if __name__=="__main__":
    name="../../binary_microchannel/Sailfish/Grid/Results/202/grid400000.npz"
    
    read(name)
