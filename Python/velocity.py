#!/usr/bin/python
import os
import numpy
import pylab
from enthought.tvtk.api import tvtk
from enthought.mayavi import mlab
import math

def Read(name):
    
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

def Read_VTI(name):
 

    gridreader = tvtk.XMLStructuredGridReader()
    gridreader.file_name = name
    gridreader.update()

    grid  = gridreader.output
    data  = grid.point_data
    points=grid.points
    dims  =grid.dimensions
    dims=dims.tolist()
    dims.reverse()
    phase= numpy.array(data.get_array("phi"))
    velocity=numpy.array(data.get_array("v"))
    velx=velocity[:,0]
    vely=velocity[:,1]
    velz=velocity[:,2]
    phase_numpy=phase.reshape(dims)
    velx_numpy =velx.reshape(dims)
    vely_numpy =vely.reshape(dims)
    velz_numpy =velz.reshape(dims)
    
    fig=mlab.figure()
    src = mlab.pipeline.scalar_field(phase_numpy)
    
    v= mlab.pipeline.vector_field(velx_numpy,vely_numpy,velz_numpy)
    vx=mlab.pipeline.scalar_field(velx_numpy)
    vy=mlab.pipeline.scalar_field(vely_numpy)
    vz=mlab.pipeline.scalar_field(velz_numpy)
    
    
    extract=mlab.pipeline.extract_grid(src)
    extract.set(z_min=1,z_max=dims[2]-2,y_min=1,y_max=dims[1]-2)
    surf = mlab.pipeline.contour_surface(extract)
 
    mlab.pipeline.image_plane_widget(vx, plane_orientation='x_axes', slice_index=250)
    #mlab.pipeline.vectors(v, mask_points=20, scale_factor=3.)
    mlab.pipeline.vector_cut_plane(v, mask_points=2, scale_factor=3)
    mlab.show()

def Read_Phase(name):
    
    print os.getcwd()
    gridreader = tvtk.XMLImageDataReader()
    gridreader.file_name = name
    gridreader.update()

    grid  = gridreader.output
    data  = grid.point_data
    dims  =grid.dimensions
    dims=dims.tolist()
    dims.reverse()
    phase= numpy.array(data.get_array("phi"))
    velocity=numpy.array(data.get_array("v"))
    velx=velocity[:,0]
    vely=velocity[:,1]
    velz=velocity[:,2]
    phase_numpy=phase.reshape(dims)
    velx_numpy =velx.reshape(dims)
    vely_numpy =vely.reshape(dims)
    velz_numpy =velz.reshape(dims)
    return phase_numpy,velx_numpy,vely_numpy,velz_numpy

    
def Show_Phase(phase_numpy,velx_numpy,vely_numpy,velz_numpy):
    dims=phase_numpy.shape
    
    fig=mlab.figure()
    src = mlab.pipeline.scalar_field(phase_numpy)
    v= mlab.pipeline.vector_field(velx_numpy,vely_numpy,velz_numpy)
    vx=mlab.pipeline.scalar_field(velx_numpy)
    vy=mlab.pipeline.scalar_field(vely_numpy)
    vz=mlab.pipeline.scalar_field(velz_numpy)
  
 
    mlab.outline()
    mlab.orientation_axes()
    extract=mlab.pipeline.extract_grid(src)
    extract.set(x_min=1,x_max=dims[0]-2,y_min=1,y_max=dims[1]-2)
    surf = mlab.pipeline.contour_surface(extract)

    v_numpy=zip(velx_numpy,vely_numpy,velz_numpy) 
    mlab.pipeline.image_plane_widget(vx, plane_orientation='x_axes', slice_index=250)
    #mlab.pipeline.vectors(v, mask_points=20, scale_factor=3.)
    mlab.pipeline.vector_cut_plane(v, mask_points=2, scale_factor=3)

    
def Analyze_Phase(name):
    
    #parameters of the binary liquid model
    k=0.04
    a=0.04
    
    phase_numpy,velx_numpy,vely_numpy,velz_numpy=Read_Phase(name)

    Show_Phase(phase_numpy,velx_numpy,vely_numpy,velz_numpy)
    dims=phase_numpy.shape

    print dims
    
    center=phase_numpy[dims[0]/2,dims[1]/2,:]
    z1 = numpy.min(numpy.where(center < 0.0))
    z2 = numpy.max(numpy.where(center < 0.0))
    if z1==0:
        z2=numpy.min(numpy.where(center>0.0))+dims[2]
        z1=numpy.max(numpy.where(center>0.0))
    pylab.figure()
    pylab.plot(center)
    print z1,z2
    
    mid =((z1+z2)/2)%dims[2]
    print mid
    slice=phase_numpy[:,:,mid]
    shape=slice.shape
    prof_axis=slice[shape[0]/2,:]
    prof_diag=numpy.diag(slice)
    
    #axis_zero=(1.0-Get_Zero(prof_axis)/(0.5*(shape[1]-2)))
    #diag_zero=(math.sqrt(2.0)-math.sqrt(2.0)*Get_Zero(prof_diag)/(0.5*(shape[1]-2)))
    
    print "Velx",velx_numpy[dims[0]/2,dims[1]/2,mid]
    print "Vely",vely_numpy[dims[0]/2,dims[1]/2,mid]
    print "Velz",velz_numpy[dims[0]/2,dims[1]/2,mid]
    
    #Calculation of the capillary number
    capillary=velx_numpy[dims[0]/2,dims[1]/2,mid]*(2.0/3.0)/math.sqrt(8.0*k*a/9.0)    
    
#    pylab.figure()
#    pylab.imshow(phase_numpy[:,:,mid])
#    pylab.figure()
#    pylab.plot(prof_axis)
#    pylab.figure()
#    pylab.plot(prof_diag)

    mlab.show()
    
    #print axis_zero,diag_zero
    print "Capillary=",capillary

    #return axis_zero,diag_zero,capillary

def Compare(name1,name2):
    import math
    arr=numpy.load(name1)

    dims=arr['phi'].shape
    velx=arr['v'][0,:,:]
    vely=arr['v'][1,:,:]
    phi=arr['phi']

    center=phi[dims[0]/2,:]
    z1 = numpy.min(numpy.where(center < 0.0))
    z2 = numpy.max(numpy.where(center < 0.0))
    if z1==0:
        z2=numpy.min(numpy.where(center>0.0))+dims[1]
        z1=numpy.max(numpy.where(center>0.0))
    pylab.figure()
    pylab.plot(center)
    print z1,z2

    mid =((z1+z2)/2)%dims[1]
    print mid
    print "Velx=",velx[dims[0]/2,mid] 
    print "Vely=",vely[dims[0]/2,mid]
    print "Dims=",dims
    print "Capillary=",velx[dims[0]/2,mid]*2.0/3.0/(math.sqrt(8.0*0.04*0.04/9.0))

    x,y=numpy.mgrid[0:dims[0],0:dims[1]]
    x_short=x[::10,::150]
    y_short=y[::10,::150]
    velx_short=velx[::10,::150]
    vely_short=vely[::10,::150]

    pylab.figure()
    pylab.imshow(arr['phi'])

    pylab.figure(figsize=[31,2.5])
    pylab.quiver(y_short,x_short,velx_short,vely_short,scale=0.4,width=0.002,headwidth=2)# ,scale=0.1)


    #3d image
    vy_3d=numpy.loadtxt(name2+"vy.txt")
    vz_3d=numpy.loadtxt(name2+"vz.txt")
    phase_3d=numpy.loadtxt(name2+"phase.txt")
    dims=phase_3d.shape
    center=phase_3d[dims[0]/2,:]
    z1 = numpy.min(numpy.where(center < 0.0))
    z2 = numpy.max(numpy.where(center < 0.0))
    if z1==0:
        z2=numpy.min(numpy.where(center>0.0))+dims[1]
        z1=numpy.max(numpy.where(center>0.0))
    pylab.figure()
    pylab.plot(center)
    print z1,z2

    mid =((z1+z2)/2)%dims[1]
    print mid
    print "Velx=",vz_3d[0,mid] 
    print "Vely=",vy_3d[0,mid]
    print "Dims=",dims
    print "Capillary=",vz_3d[0,mid]*2.0/3.0/(math.sqrt(8.0*0.04*0.04/9.0))

    x,y=numpy.mgrid[0:dims[0],0:dims[1]]
    x_short=x[::3,::75]
    y_short=y[::3,::75]
    velx_short=vz_3d[::3,::75]
    vely_short=vy_3d[::3,::75]

    pylab.figure()
    pylab.imshow(phase_3d, origin="lower")

    pylab.figure(figsize=[31,2.5])
    pylab.quiver(y_short,x_short,velx_short-0.013,vely_short,scale=0.6,width=0.0023,headwidth=2)# ,scale=0.1)
    pylab.contour(phase_3d,[0.0])
    pylab.show()


if __name__=="__main__":
    #name1="../../binary_microchannel/Sailfish/Grid/Results/202/grid400000.npz"
    #Mac filenames
    #name1="20/capillary200000.npz"
    #name2="4/"
    
    #Linux filenames
    name1="../../binary_microchannel/Sailfish/Capillary/Results/20/capillary200000.npz"
    name2="../Results/Force0000002/4/"
    #Analyze_Phase(name)
    Compare(name1,name2)
