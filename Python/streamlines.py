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

def read_vtk_2d(name):
    import vtk
    import numpy

    gridreader = vtk.vtkXMLStructuredGridReader()
    gridreader.SetFileName(name)
    #gridreader.SetPointArrayStatus("Density",0)
    selection=gridreader.GetPointDataArraySelection()
    selection.DisableArray("Density")
    #selection.DisableArray("Velocity")
    #selection.DisableArray("Phase")
    gridreader.Update()
    
    grid  = gridreader.GetOutput()

    data  = grid.GetPointData()
    points=grid.GetPoints()
    dims  =grid.GetDimensions()
    #data.SetActiveScalars("Phase"); 
    data.SetActiveVectors("Velocity")
    velocity=data.GetArray("Velocity")
    phase = data.GetArray("Phase")
    
    vel_image=vtk.vtkImageData()
    vel_image.SetDimensions(dims[1],dims[2],1)
    vel_image.SetSpacing(1.0, 1.0, 1.0)
    vel_image.SetOrigin(0.0, 0.0, 0.0)

    
    Array = vtk.vtkDoubleArray()
    Array.SetNumberOfComponents(3)
    Array.SetNumberOfTuples(points.GetNumberOfPoints())
    Array.Initialize()
    Array.SetName("VelocitySlice")

    ArrayPhase = vtk.vtkDoubleArray()
    ArrayPhase.SetNumberOfComponents(1)
    ArrayPhase.SetNumberOfTuples(points.GetNumberOfPoints())
    ArrayPhase.Initialize()
    ArrayPhase.SetName("PhaseSlice")
    
    for counter in range(0,10):
        print points.GetPoint(counter)
    raw_input()
    exam=1
    print dims[0],dims[1],dims[2]
    #raw_input()
    for counter in range(0,vel_image.GetNumberOfPoints()):
        x,y,z=vel_image.GetPoint(counter)
        counter_big=y*dims[0]*dims[1]+x*dims[0]+exam
        #print x,y,z
        #print counter_big
        
        vz=velocity.GetTuple3(counter_big)[2]
        vy=velocity.GetTuple3(counter_big)[1]
        phase_value=phase.GetTuple1(counter_big)
        Array.InsertNextTuple3(vy,vz,0.0)
        ArrayPhase.InsertNextTuple1(phase_value)
        #if phase_value<0.0:
        #    print phase_value
    
    vel_image.GetPointData().SetVectors(Array)
    vel_image.GetPointData().SetScalars(ArrayPhase)
    vel_image.GetPointData().SetActiveScalars("PhaseSlice")
    vel_image.GetPointData().SetActiveVectors("VelocitySlice")
    print vel_image.GetPointData().GetArray("PhaseSlice")
    vel_image.Update()
    
    
    
    contour=vtk.vtkContourFilter()
    contour.SetInput(vel_image)
    #contour.SetValue(0,0.0)

    contourMapper = vtk.vtkPolyDataMapper()
    #contourMapper.SetScalarRange(phase.GetRange())
    contourMapper.SetInputConnection(contour.GetOutputPort())
    #contourMapper.SetColorMode
 
    contourActor = vtk.vtkActor()
    contourActor.SetMapper(contourMapper)

    #contourActor.SetColor(0.0,0.1,0.2)
    
    rk4=vtk.vtkRungeKutta4()
    line = vtk.vtkLineSource()
    line.SetPoint1(50, 0, 0)
    line.SetPoint2(0, 0, 0)
    line.SetResolution(21)
    streamer = vtk.vtkStreamLine()
    streamer.SetInput(vel_image)
    streamer.SetSource(line.GetOutput())
    streamer.SetMaximumPropagationTime(300000)
    #streamer.GetStepLengthMinValue()=0.01
    #streamer.GetStepLengthMaxValue()=0.5
    #streamer.SetTerminalSpeed(1e-13)
    streamer.SetIntegrationStepLength(.2)
    streamer.SetStepLength(.5)
    #streamer.SetNumberOfThreads(1)
    streamer.SetIntegrationDirectionToIntegrateBothDirections() 
    #streamer.VorticityOn()
    streamer.SetIntegrator(rk4)

    streamMapper = vtk.vtkPolyDataMapper()
    streamMapper.SetInputConnection(streamer.GetOutputPort())
    #streamMapper.SetScalarRange(data.GetOutput().GetScalarRange())
    streamActor = vtk.vtkActor()
    streamActor.SetMapper(streamMapper)
    streamActor.VisibilityOn()

    
    # Create the Renderer, RenderWindow, and RenderWindowInteractor
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    # Add the actors to the render; set the background and size
    ren.AddActor(streamActor)
    ren.AddActor(contourActor)
    #ren.SetBackground(1.0,1.0,1.0)
    ren.SetBackground(0.1, 0.2, 0.4)
    
    renWin.SetSize(500, 500)

    #ren.ResetCamera()
    #cam1 = ren.GetActiveCamera()
    #cam1.SetPosition(1000,0,500)
    #cam1.SetFocalPoint(0,0,500)

    iren.Initialize()
    renWin.Render()
    iren.Start()


    
def read_vtk(name):
    import vtk
    import numpy

    gridreader = vtk.vtkXMLStructuredGridReader()
    gridreader.SetFileName(name)
    #gridreader.SetPointArrayStatus("Density",0)
    selection=gridreader.GetPointDataArraySelection()
    selection.DisableArray("Density")
    #selection.DisableArray("Velocity")
    selection.DisableArray("Phase")
    gridreader.Update()
    
    grid  = gridreader.GetOutput()

    data  = grid.GetPointData()
    points=grid.GetPoints()
    dims  =grid.GetDimensions()
    #data.SetActiveScalars("Phase"); 
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



    plane = vtk.vtkPlane()
    #plane.SetInput(grid)
    #plane.SetResolution(50, 50)
    plane.SetOrigin(25.0,10.0,0.0)
    plane.SetNormal(1.0,0.0,0.0)

    cutter=vtk.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInput(grid)
    cutter.Update()
 
    cutterMapper =vtk.vtkPolyDataMapper()
    cutterMapper.SetInputConnection(cutter.GetOutputPort(0))
 
    # Create plane actor
    planeActor =vtk.vtkActor()
    planeActor.GetProperty().SetColor(1.0,1,0);
    planeActor.GetProperty().SetLineWidth(2);
    planeActor.SetMapper(cutterMapper);

    #transP1 = vtk.vtkTransform()
    #transP1.Translate(3.7, 0.0, 28.37)
    #transP1.Scale(5, 5, 5)
    #transP1.RotateY(90)
    #tpd1 = vtk.vtkTransformPolyDataFilter()
    #tpd1.SetInputConnection(plane.GetOutputPort())
    #tpd1.SetTransform(transP1)
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
    #ren.AddActor(contourActor)
    #ren.AddActor(outlineActor)
    ren.AddActor(streamActor)
    #ren.AddActor(rake1Actor)
    #ren.AddActor(rake2Actor)
    #ren.AddActor(planeActor)h=streamline(vz,vy,sz,sy,[0.1, 15000])
    #ren.SetBackground(1.0,1.0,1.0)
    ren.SetBackground(0.1, 0.2, 0.4)
    
    renWin.SetSize(500, 500)

    # Zoom in closer
    #ren.ResetCamera()
    cam1 = ren.GetActiveCamera()
    cam1.SetPosition(1000,0,500)
    cam1.SetFocalPoint(0,0,500)

    iren.Initialize()
    renWin.Render()
    iren.Start()
def draw_streamlines(name):
    import vtk
    import numpy
    from PyNGL import Ngl

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
    velocity=data.GetArray("Velocity")
    phase=data.GetArray("Phase")
    vz=numpy.zeros([dims[1],dims[2]])
    vy=numpy.zeros([dims[1],dims[2]])
    phase_numpy=numpy.zeros([dims[1],dims[2]])
    print vz.shape
    print vy.shape
    for counter in range(0,points.GetNumberOfPoints()):
        coorx,coory,coorz=points.GetPoint(counter)
        #print coorx,coory,coorz
        if coorx==1:
            vz[int(coory),int(coorz)]=velocity.GetTuple3(counter)[2]
            vy[int(coory),int(coorz)]=velocity.GetTuple3(counter)[1]
            phase_numpy[int(coory),int(coorz)]=phase.GetTuple1(counter)
    #data.SetActiveScalars("Phase"); 
    #data.SetActiveVectors("Velocity")
    numpy.savetxt("vz.txt",vz)
    numpy.savetxt("vy.txt",vy)
    numpy.savetxt("phase.txt",phase_numpy)
    wks_type = "ps"
    wks = Ngl.open_wks(wks_type,"test")
    resources = Ngl.Resources()
     
     #uvar = file.variables["U_GRD_6_ISBL"]
     #vvar = file.variables["V_GRD_6_ISBL"]
     #if hasattr(uvar,"units"):
     #  resources.tiMainString = "GRD_6_ISBL (u,v " + uvar.units + ")"
     #else:
        #resources.tiMainString = "GRD_6_ISBL"
     #if hasattr(uvar,"_FillValue"):
     #    resources.vfMissingUValueV = uvar._FillValue
     # if hasattr(vvar,"_FillValue"):
     #    resources.vfMissingVValueV = vvar._FillValue
  
    resources.tiMainFont    = "Times-Roman"
    resources.tiXAxisString = "streamlines"
    resources.vpHeightF = 0.25 # Define height, width, and location of plot.
    resources.vpWidthF  = 7.5*0.25

     
     #plot = Ngl.streamline(wks,uvar[0,::2,::2],vvar[0,::2,::2],resources) 
    plot=Ngl.streamline(wks,vz[:,::10],vy[:,::10],resources)
     
    #Ngl.end()

    
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
    #read_vtk(name)
    read_vtk_2d(name)
    #draw_streamlines(name)
    #show(name,300)