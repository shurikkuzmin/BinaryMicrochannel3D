#!/usr/bin/python
def read(name):
    from enthought.tvtk.api import tvtk
    from enthought.mayavi import mlab
    import numpy


    gridreader = tvtk.XMLStructuredGridReader()
    gridreader.file_name = name
    gridreader.update()

    grid  = gridreader.output
    data  = grid.point_data
    points=grid.points
    dims  =grid.dimensions
    dims=dims.tolist()
    dims.reverse()
    phase= numpy.array(data.get_array("Phase"))
    velocity=numpy.array(data.get_array("Velocity"))
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
 
    mlab.axes()
    mlab.pipeline.image_plane_widget(vx, plane_orientation='x_axes', slice_index=250)
    #mlab.pipeline.vectors(v, mask_points=20, scale_factor=3.)
    mlab.pipeline.vector_cut_plane(v, mask_points=2, scale_factor=3)
    mlab.show()

def show(name,slice):
    from enthought.tvtk.api import tvtk
    import pylab
    import numpy


    gridreader = tvtk.XMLStructuredGridReader()
    gridreader.file_name = name
    gridreader.update()

    grid  = gridreader.output
    data  = grid.point_data
    points=grid.points
    dims  =grid.dimensions
    dims=dims.tolist()
    dims.reverse()
    phase= numpy.array(data.get_array("Phase"))
    velocity=numpy.array(data.get_array("Velocity"))
    velx=velocity[:,0]
    vely=velocity[:,1]
    velz=velocity[:,2]
    phase_numpy=phase.reshape(dims)
    velx_numpy =velx.reshape(dims)
    vely_numpy =vely.reshape(dims)
    velz_numpy =velz.reshape(dims)
    
    pylab.imshow(phase_numpy[slice,:,:])
    pylab.show()


def check_symmetry(name):
    from enthought.tvtk.api import tvtk
    import pylab
    import numpy


    gridreader = tvtk.XMLStructuredGridReader()
    gridreader.file_name = name
    gridreader.update()

    grid  = gridreader.output
    data  = grid.point_data
    points=grid.points
    dims  =grid.dimensions
    dims=dims.tolist()
    dims.reverse()
    phase= numpy.array(data.get_array("Phase"))
    velocity=numpy.array(data.get_array("Velocity"))
    velx=velocity[:,0]
    vely=velocity[:,1]
    velz=velocity[:,2]
    phase_numpy=phase.reshape(dims)
    velx_numpy =velx.reshape(dims)
    vely_numpy =vely.reshape(dims)
    velz_numpy =velz.reshape(dims)

    print velx_numpy[25,:10,0], velx_numpy[25,:10,1]
    print vely_numpy[25,:10,0], vely_numpy[25,:10,1]
    print velz_numpy[25,:10,0], velz_numpy[25,:10,1]

    print "Another set:"
    print velx_numpy[25,0,:10], velx_numpy[25,1,:10]
    print vely_numpy[25,0,:10], vely_numpy[25,1,:10]
    print velz_numpy[25,0,:10], velz_numpy[25,1,:10]
    
    pylab.figure()
    pylab.imshow(phase_numpy[25,:,:])
    pylab.figure()
    pylab.imshow(velz_numpy[25,:,:])
    pylab.show()
    print velx_numpy.shape


def compare(name,name_quarter):
    from enthought.tvtk.api import tvtk
    import pylab
    import numpy
    
    gridreader = tvtk.XMLStructuredGridReader()
    gridreader.file_name = name
    gridreader.update()

    gridreader_quarter = tvtk.XMLStructuredGridReader()
    gridreader_quarter.file_name = name_quarter
    gridreader_quarter.update()

    grid  = gridreader.output
    data  = grid.point_data
    
    grid_quarter=gridreader_quarter.output
    data_quarter = grid_quarter.point_data
    
    dims  =grid.dimensions
    dims=dims.tolist()
    dims.reverse()

    dims_quarter  =grid_quarter.dimensions
    dims_quarter=dims_quarter.tolist()
    dims_quarter.reverse()
    
    phase= numpy.array(data.get_array("Phase"))
    velocity=numpy.array(data.get_array("Velocity"))
    velx=velocity[:,0]
    vely=velocity[:,1]
    velz=velocity[:,2]
    phase_numpy=phase.reshape(dims)
    velx_numpy =velx.reshape(dims)
    vely_numpy =vely.reshape(dims)
    velz_numpy =velz.reshape(dims)
 
    phase_quarter= numpy.array(data_quarter.get_array("Phase"))
    velocity_quarter=numpy.array(data_quarter.get_array("Velocity"))
    velx_quarter=velocity_quarter[:,0]
    vely_quarter=velocity_quarter[:,1]
    velz_quarter=velocity_quarter[:,2]
    phase_quarter_numpy=phase_quarter.reshape(dims_quarter)
    velx_quarter_numpy =velx_quarter.reshape(dims_quarter)
    vely_quarter_numpy =vely_quarter.reshape(dims_quarter)
    velz_quarter_numpy =velz_quarter.reshape(dims_quarter)

    pylab.figure()
    pylab.imshow(phase_numpy[25,0:26,0:26])
    pylab.figure()
    pylab.imshow(phase_quarter_numpy[25,:,:],origin="lower")
    pylab.figure()
    pylab.plot(phase_quarter_numpy[25,:,0])
    pylab.plot(phase_quarter_numpy[25,:,1])
    
    print velx_quarter_numpy[25,:10,0], velx_quarter_numpy[25,:10,1]
    print vely_quarter_numpy[25,:10,0], vely_quarter_numpy[25,:10,1]
    print velz_quarter_numpy[25,:10,0], velz_quarter_numpy[25,:10,1]
    
    pylab.show()
 
 
    

if __name__=="__main__":
    name_symmetry="../QuarterCode/phase00100.vts"
    name="../Code/phase00300.vts"
    name_quarter="../QuarterCode/phase00300.vts"
    #read(name)
    #show(n,300)
    #compare(name,name_quarter)
    check_symmetry(name_symmetry)