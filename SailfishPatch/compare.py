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

if __name__=="__main__":
    name="test1000_0.vti"
    read_michal(name)
    #show(name,300)