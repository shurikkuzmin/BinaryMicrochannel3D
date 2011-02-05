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
def read_paraview(name):
    #import paraview.simple
    #from paraview.simple import *
    
    reader=XMLStructuredGridReader(FileName=name)
    reader.UpdatePipeline()
    
    contour=Contour()
    contour.ContourBy="Phase"
    contour.Isosurfaces=[0.0]
    Show(contour)
    Render()
    #Show()
    Start()
    
def read_vtk(name):
    import vtk
    import numpy

    gridreader = vtk.vtkXMLStructuredGridReader()
    gridreader.SetFileName(name)
    #gridreader.SetPointArrayStatus("Density",0)
    #selection=gridreader.GetPointDataArraySelection()
    #selection.DisableArray("Density")
    #selection.DisableArray("Velocity")
    gridreader.Update()
    
    grid  = gridreader.GetOutput()

    data  = grid.GetPointData()
    points=grid.GetPoints()
    dims  =grid.GetDimensions()
    data.SetActiveScalars("Phase"); 
    data.SetActiveVectors("Velocity")
    #Create a contour    
    contour=vtk.vtkContourFilter()
    contour.SetInputConnection(gridreader.GetOutputPort())
    contour.SetValue(0,0.0)

    contourMapper = vtk.vtkPolyDataMapper()
    #contourMapper.SetScalarRange(phase.GetRange())
    contourMapper.SetInputConnection(contour.GetOutputPort())
 
    contourActor = vtk.vtkActor()
    contourActor.SetMapper(contourMapper)

    rake1 = vtk.vtkLineSource()
    rake1.SetPoint1(0, 0, 0)
    rake1.SetPoint2(0, 50, 0)
    rake1.SetResolution(21)
    rake1Mapper = vtk.vtkPolyDataMapper()
    rake1Mapper.SetInputConnection(rake1.GetOutputPort())
    rake1Actor = vtk.vtkActor()
    rake1Actor.SetMapper(rake1Mapper)

    rake2 = vtk.vtkLineSource()
    rake2.SetPoint1(0, 0, 850)
    rake2.SetPoint2(0, 50,850)
    rake2.SetResolution(21)
    rake2Mapper = vtk.vtkPolyDataMapper()
    rake2Mapper.SetInputConnection(rake2.GetOutputPort())
    rake2Actor = vtk.vtkActor()
    rake2Actor.SetMapper(rake2Mapper)

    lines=vtk.vtkAppendPolyData()
    lines.AddInput(rake1.GetOutput())
    lines.AddInput(rake2.GetOutput())
    rk4 = vtk.vtkRungeKutta4()
    streamer = vtk.vtkStreamLine()
    streamer.SetInputConnection(gridreader.GetOutputPort())
    streamer.SetSource(lines.GetOutput())
    streamer.SetMaximumPropagationTime(100000)
    #streamer.GetStepLengthMinValue()=0.01
    #streamer.GetStepLengthMaxValue()=0.5
    #streamer.SetTerminalSpeed(1e-13)
    streamer.SetIntegrationStepLength(.2)
    streamer.SetStepLength(.2)
    #streamer.SetNumberOfThreads(1)
    streamer.SetIntegrationDirectionToIntegrateBothDirections() 
    #streamer.VorticityOn()
    streamer.SetIntegrator(rk4)

    #integ = vtk.vtkRungeKutta4()
    #sl = vtk.vtkStreamLine()
    #sl.SetInputConnection(pl3d.GetOutputPort())
    #sl.SetSource(rake.GetOutput())
    #sl.SetIntegrator(integ)
    #sl.SetMaximumPropagationTime(0.1)
    #sl.SetIntegrationStepLength(0.1)
    #sl.SetIntegrationDirectionToBackward()
    #sl.SetStepLength(0.001)

    #plane = vtk.vtkPlaneSource()
    #plane.SetResolution(50, 50)
    #transP1 = vtk.vtkTransform()
    #transP1.Translate(3.7, 0.0, 28.37)
    #transP1.Scale(5, 5, 5)
    #transP1.RotateY(90)
    #tpd1 = vtk.vtkTransformPolyDataFilter()
    #tpd1.SetInputConnection(plane.GetOutputPort())
    #tpd1.SetTransform(transP1)
    #outTpd1 = vtk.vtkOutlineFilter()
    #outTpd1.SetInputConnection(tpd1.GetOutputPort())
    #mapTpd1 = vtk.vtkPolyDataMapper()
    #mapTpd1.SetInputConnection(outTpd1.GetOutputPort())
    #tpd1Actor = vtk.vtkActor()
    #tpd1Actor.SetMapper(mapTpd1)
    #tpd1Actor.GetProperty().SetColor(0, 0, 0)

    #lineWidget = vtk.vtkLineWidget()
    #seeds = vtk.vtkPolyData()
    #lineWidget.SetInput(grid)
    #lineWidget.SetAlignToYAxis()
    #lineWidget.PlaceWidget()
    #lineWidget.GetPolyData(seeds)
    #lineWidget.ClampToBoundsOn()

    
    streamMapper = vtk.vtkPolyDataMapper()
    streamMapper.SetInputConnection(streamer.GetOutputPort())
    #streamMapper.SetScalarRange(data.GetOutput().GetScalarRange())
    streamActor = vtk.vtkActor()
    streamActor.SetMapper(streamMapper)
    streamActor.VisibilityOn()


    #Create the outline
    outline = vtk.vtkOutlineFilter()
    outline.SetInputConnection(gridreader.GetOutputPort())

    outlineMapper = vtk.vtkPolyDataMapper()
    outlineMapper.SetInputConnection(outline.GetOutputPort())

    outlineActor = vtk.vtkActor()
    outlineActor.SetMapper(outlineMapper)
    #outlineActor.GetProperty().SetColor(0,0,0)

   
    # Create the Renderer, RenderWindow, and RenderWindowInteractor
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    # Add the actors to the render; set the background and size
    ren.AddActor(contourActor)
    ren.AddActor(outlineActor)
    ren.AddActor(streamActor)
    ren.AddActor(rake1Actor)
    ren.AddActor(rake2Actor)
    #ren.SetBackground(1.0,1.0,1.0)
    ren.SetBackground(0.1, 0.2, 0.4)
    renWin.SetSize(500, 500)

    # Zoom in closer
    ren.ResetCamera()
    cam1 = ren.GetActiveCamera()
    cam1.Zoom(1.4)

    iren.Initialize()
    renWin.Render()
    iren.Start()


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
    #name="/Users/shurik/Documents/Temp/Paraview/phase10000.vts"
    #read_paraview(name)
    read_vtk(name)
    #show(name,300)