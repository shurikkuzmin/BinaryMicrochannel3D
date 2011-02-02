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
    print velocity.shape
    velx=velocity[:,0]
    vely=velocity[:,1]
    velz=velocity[:,2]
    phase_numpy=phase.reshape(dims)
    velx_numpy =velx.reshape(dims)
    vely_numpy =vely.reshape(dims)
    velz_numpy =velz.reshape(dims)
    
    fig=mlab.figure()
    src = mlab.pipeline.scalar_field(phase_numpy)
    v= mlab.pipeline.vector_field(velz_numpy,vely_numpy,velx_numpy)
    #vx=mlab.pipeline.scalar_field(velx_numpy)
    #vy=mlab.pipeline.scalar_field(vely_numpy)
    #vz=mlab.pipeline.scalar_field(velz_numpy)
    
    #extract=mlab.pipeline.extract_grid(src)
    #extract.set(z_min=1,z_max=dims[2]-2,y_min=1,y_max=dims[1]-2)
    #surf = mlab.pipeline.contour_surface(extract)
    surf = mlab.pipeline.contour_surface(src)
 
    #mlab.pipeline.image_plane_widget(vx, plane_orientation='x_axes', slice_index=250)
    #mlab.pipeline.vectors(v, mask_points=20, scale_factor=3.)
    mlab.pipeline.vector_cut_plane(v,plane_orientation='y_axes', mask_points=5)
    mlab.show()
def read_vtk(name):
    import vtk
    import numpy

    gridreader = vtk.vtkXMLStructuredGridReader()
    gridreader.SetFileName(name)
    gridreader.Update()

    grid  = gridreader.GetOutput()
    data  = grid.GetPointData()
    points=grid.GetPoints()
    dims  =grid.GetDimensions()
    
    
    #print dims
    #dims=list(dims)
    #dims.reverse()
    

    phase= data.GetArray("Phase")
    #velocity=data.GetArray("Velocity")
    
    contour=vtk.vtkContourFilter()
    contour.SetInputConnection(gridreader.GetOutputPort())
    contour.GenerateValues(50, phase.GetRange())
    
    contourMapper = vtk.vtkPolyDataMapper()
    contourMapper.SetInputConnection(contour.GetOutputPort())
    #contourMapper.SetScalarRange(gridreader.GetOutput().GetScalarRange())

    #sr=vtk.vtkStructuredGridToPolyDataFilter(gridreader.GetOutputPort())
    stlMapper = vtk.vtkPolyDataMapper()
    #stlMapper.SetInput(data)
    stlMapper.SetInputConnection(contour.GetOutputPort())

    contour=vtk.vtkContourFilter()
    
    #reslice.SetInputConnection(reader.GetOutputPort())
    #reslice.SetOutputDimensionality(2)
    #reslice.SetResliceAxes(sagittal)
    #reslice.SetInterpolationModeToLinear()



    stlActor = vtk.vtkActor()
    stlActor.SetMapper(stlMapper)

    # Create the Renderer, RenderWindow, and RenderWindowInteractor
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    # Add the actors to the render; set the background and size
    ren.AddActor(stlActor)
    ren.SetBackground(0.1, 0.2, 0.4)
    renWin.SetSize(500, 500)

    # Zoom in closer
    ren.ResetCamera()
    cam1 = ren.GetActiveCamera()
    cam1.Zoom(1.4)

    iren.Initialize()
    renWin.Render()
    iren.Start()

    
    #phase_array=[]
    #velocity_array=[]
    #phase.GetData(0,phase.GetNumberOfTuples(),0,0,phase_array)
    #velocity.GetData(0,velocity.GetNumberOfTuples(),0,0,velocity_array)
    
    
    
    
    #print velocity.shape
    #print data.GetArrayName(0)
    #print data.GetArrayName(1)
    #print data.GetArrayName(2)
    
    #print velocity
    #velx=velocity[:,0]
    #vely=velocity[:,1]
    #velz=velocity[:,2]
    
    #phase_numpy=phase.reshape(dims)
    #velx_numpy =velx.reshape(dims)
    #vely_numpy =vely.reshape(dims)
    #velz_numpy =velz.reshape(dims)
    
    
    

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
    name="../Temp/phase040000.vts"
    read_vtk(name)
    #show(name,300)