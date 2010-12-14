#!/usr/bin/python
def read_my(name):
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
 
    mlab.pipeline.image_plane_widget(vx, plane_orientation='x_axes', slice_index=250)
    #mlab.pipeline.vectors(v, mask_points=20, scale_factor=3.)
    mlab.pipeline.vector_cut_plane(v, mask_points=2, scale_factor=3)
    mlab.show()

def read_michal(name):
    from enthought.tvtk.api import tvtk
    from enthought.mayavi import mlab
    import numpy


    gridreader = tvtk.XMLImageDataReader()
    gridreader.file_name = name
    gridreader.update()

    grid  = gridreader.output
    data  = grid.point_data
    #points=grid.points
    dims  =grid.dimensions
    dims=dims.tolist()
    dims.reverse()
    phase= numpy.array(data.get_array("phi"))
    #velocity=numpy.array(data.get_array("Velocity"))
    #velx=velocity[:,0]
    #vely=velocity[:,1]
    #velz=velocity[:,2]
    phase_numpy=phase.reshape(dims)
    #velx_numpy =velx.reshape(dims)
    #vely_numpy =vely.reshape(dims)
    #velz_numpy =velz.reshape(dims)
    
    fig=mlab.figure()
    src = mlab.pipeline.scalar_field(phase_numpy)
    #v= mlab.pipeline.vector_field(velx_numpy,vely_numpy,velz_numpy)
    #vx=mlab.pipeline.scalar_field(velx_numpy)
    #vy=mlab.pipeline.scalar_field(vely_numpy)
    #vz=mlab.pipeline.scalar_field(velz_numpy)
    
    extract=mlab.pipeline.extract_grid(src)
    extract.set(z_min=1,z_max=dims[2]-2,y_min=1,y_max=dims[1]-2)
    surf = mlab.pipeline.contour_surface(extract)
 
    #mlab.pipeline.image_plane_widget(vx, plane_orientation='x_axes', slice_index=250)
    #mlab.pipeline.vectors(v, mask_points=20, scale_factor=3.)
    #mlab.pipeline.vector_cut_plane(v, mask_points=2, scale_factor=3)
    mlab.show()

def get_zero(prof):
    zero=0
    #pylab.figure()
    #pylab.plot(prof)
    for counter in range(0, len(prof)/2):
        if prof[counter]>=0 and prof[counter+1]<0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    return zero-0.5

def show(name,slice):
    from enthought.tvtk.api import tvtk
    import pylab
    import numpy


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
    arr=phase_numpy[:,:,slice]
    
    pylab.figure()
    pylab.imshow(arr)
    pylab.figure()
    pylab.plot(numpy.diag(arr))
    print get_zero(numpy.diag(arr))
    print velz_numpy[27,27,slice]
    print velx_numpy[27,27,slice]
    print vely_numpy[27,27,slice]
    pylab.show()

def compare(name_michal,name_alex,slice):
    from enthought.tvtk.api import tvtk
    import pylab
    import numpy

    gridreader_alex = tvtk.XMLStructuredGridReader()
    gridreader_alex.file_name = name_alex
    gridreader_alex.update()

    gridreader_michal = tvtk.XMLImageDataReader()
    gridreader_michal.file_name = name_michal
    gridreader_michal.update()

    grid_alex  = gridreader_alex.output
    data_alex  = grid_alex.point_data
    dims_alex  = grid_alex.dimensions
    dims_alex  = dims_alex.tolist()
    dims_alex.reverse()
    phase_alex = numpy.array(data_alex.get_array("Phase"))
    phase_alex = phase_alex.reshape(dims_alex)

    grid_michal  = gridreader_michal.output
    data_michal  = grid_michal.point_data
    dims_michal  = grid_michal.dimensions
    dims_michal  = dims_michal.tolist()
    dims_michal.reverse()
    phase_michal = numpy.array(data_michal.get_array("phi"))
    phase_michal = phase_michal.reshape(dims_michal)
    
    print phase_alex.shape
    print phase_michal.shape
    
    surf_alex = phase_alex[slice,:,:]
    surf_michal=phase_michal[:,:,slice]
    pylab.figure()
    pylab.imshow(surf_alex)
    pylab.figure()
    pylab.imshow(surf_michal)
    
    pylab.figure()

    prof_michal=surf_michal[:,surf_michal.shape[1]/2]
    prof_alex=surf_alex[:,surf_alex.shape[1]/2]
    pylab.plot(prof_michal[1:-1])
    pylab.plot(prof_alex[1:-1])
    print numpy.max(numpy.abs(prof_michal-prof_alex))
    print prof_michal-prof_alex
    pylab.show()


 
if __name__=="__main__":
    name_michal="../Temp/asym100000_0.vti"
    #name_alex="/phase10000.vts"
    slice=350
    show(name_michal,slice)
    #read_my(name_alex)
    #compare(name_michal,name_alex,slice)
    #show(name,300)