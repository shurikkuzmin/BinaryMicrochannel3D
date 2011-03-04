#! /usr/bin/python
import vtk
import numpy
import math
def contour_movie(name):
    gridreader = vtk.vtkXMLStructuredGridReader()
    #gridreader=vtk.vtkImageReader()    
    gridreader.SetFileName(name)
    gridreader.Update()
    
    grid  = gridreader.GetOutput()
    data  = grid.GetPointData()
    points=grid.GetPoints()
    dims  =grid.GetDimensions()
    #velocity=data.GetArray("v")
    #phase=data.GetArray("phi")
    #print phase
    velocity=data.GetArray("Velocity")
    phase=data.GetArray("Phase")
    #data.SetActiveScalars("Phase")
    
    phase_image=vtk.vtkImageData()
    phase_image.SetSpacing(1.0,1.0,1.0)
    phase_image.SetOrigin(0.0,0.0,0.0)
    phase_image.SetDimensions(dims[0],dims[1],dims[2])
    phase_image.GetPointData().SetScalars(phase)
    phase_image.GetPointData().SetVectors(velocity)
    
    outline=vtk.vtkOutlineFilter()
    outline.SetInput(phase_image)
    
    outlineMapper=vtk.vtkDataSetMapper()
    outlineMapper.SetInputConnection(outline.GetOutputPort())
    
    outlineActor=vtk.vtkActor()
    outlineActor.SetMapper(outlineMapper)
    
    extract=vtk.vtkExtractVOI()
    extract.SetInput(phase_image)
    extract.SetVOI(0,0,0,dims[1]-1,0,dims[2]-1)
    extract.Update()
     
    contour=vtk.vtkContourFilter()
    #contour.SetInputConnection(phase_image.GetOutputPort())
    contour.SetInput(phase_image)
    contour.SetValue(0,0.0)
    contour.Update()

    contourMapper = vtk.vtkPolyDataMapper()
    #contourMapper.SetScalarRange(phase.GetRange())
    contourMapper.SetInputConnection(contour.GetOutputPort())
    
    contourActor = vtk.vtkActor()
    contourActor.SetMapper(contourMapper)
     
    ren = vtk.vtkRenderer()
    ren.AddActor(contourActor)
    ren.AddActor(outlineActor)
    ren.SetBackground(0.1, 0.2, 0.4)
    #ren.SetColor(0.1,0.5,0.2)    
    #ren.SetViewport(0, 0, .3, 1)

    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    renWin.SetSize(500, 500)
    
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
 
    iren.Initialize()
    renWin.Render()
 
    numbers=720
    delta=2.0*math.pi/numbers
    
    cam=ren.GetActiveCamera()
    #cam.SetPosition(300,0,750)
    renWin.Render()
    xinit,yinit,zinit=cam.GetPosition()
    radius=math.sqrt(xinit**2+yinit**2+zinit**2)    
    print xinit,yinit,zinit
    cam.SetViewUp(1.00000,0.0000001,0.00000)    
    cam.Zoom(3)    
    for counter,alpha in enumerate(numpy.arange(0,2.0*math.pi,delta)):
        cam.Azimuth(delta*180.0/math.pi)
        renWin.Render()
        im_filter=vtk.vtkWindowToImageFilter()
        im_filter.SetInput(renWin)
        writer=vtk.vtkJPEGWriter()
        file_name="tmp/contour"+str(0)*(3-len(str(counter)))+str(counter)+".jpg"
        writer.SetFileName(file_name)
        writer.SetInput(im_filter.GetOutput())    
        writer.Write()

    
    iren.Start()
       
def contour_movie_vti(name):
    #gridreader = vtk.vtkXMLStructuredGridReader()
    gridreader=vtk.vtkXMLImageDataReader()    
    gridreader.SetFileName(name)
    gridreader.Update()
    
    
    grid  = gridreader.GetOutput()
    data  = grid.GetPointData()
    #dims  =grid.GetDimensions()
    #print "Dims=",dims    
    velocity=data.GetArray("v")
    phase=data.GetArray("phi")
    data.SetScalars(phase)
    data.SetVectors(velocity)
    print "Phase=",phase
    print "Data=",data
    
    outline=vtk.vtkOutlineFilter()
    outline.SetInputConnection(gridreader.GetOutputPort())
    #outline.SetInput(grid)    
    
    outlineMapper=vtk.vtkDataSetMapper()
    outlineMapper.SetInputConnection(outline.GetOutputPort())
    
    outlineActor=vtk.vtkActor()
    outlineActor.SetMapper(outlineMapper)
    
     
    contour=vtk.vtkContourFilter()
    #contour.SetInputConnection(phase_image.GetOutputPort())
    #contour.SetInput(phase_image)
    contour.SetInputConnection(gridreader.GetOutputPort())    
    contour.SetValue(0,0.0)
    contour.Update()

    contourMapper = vtk.vtkPolyDataMapper()
    #contourMapper.SetScalarRange(phase.GetRange())
    contourMapper.SetInputConnection(contour.GetOutputPort())
    
    contourActor = vtk.vtkActor()
    contourActor.SetMapper(contourMapper)
     
    ren = vtk.vtkRenderer()
    ren.AddActor(contourActor)
    ren.AddActor(outlineActor)
    ren.SetBackground(0.1, 0.2, 0.4)
    #ren.SetColor(0.1,0.5,0.2)    
    #ren.SetViewport(0, 0, .3, 1)

    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    renWin.SetSize(500, 500)
    
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
 
    iren.Initialize()
    renWin.Render()
 
    numbers=720
    delta=2.0*math.pi/numbers
    
    cam=ren.GetActiveCamera()
    print "Camera=", cam.GetViewUp()
    xinit,yinit,zinit=cam.GetPosition()
    radius=math.sqrt(xinit**2+yinit**2+zinit**2)    
    print xinit,yinit,zinit
    #cam.Azimuth(delta*50)
    #contourActor.RotateWXYZ(25,1,0,0)
    #renWin.Render()
    #cam.Zoom(3)
    cam.SetViewUp(1,0,0)
    for counter,alpha in enumerate(numpy.arange(0,2.0*math.pi,delta)):
        cam.Azimuth(delta*180.0/math.pi)
        #contourActor.RotateWXYZ(delta*180/math.pi,0,1,0)
        #outlineActor.RotateWXYZ(delta*180/math.pi,0,1,0)

        renWin.Render()
        im_filter=vtk.vtkWindowToImageFilter()
        im_filter.SetInput(renWin)
        writer=vtk.vtkJPEGWriter()
        file_name="tmp/contour"+str(0)*(3-len(str(counter)))+str(counter)+".jpg"
        writer.SetFileName(file_name)
        writer.SetInput(im_filter.GetOutput())    
        writer.Write()

    
    iren.Start()
    
   
    
    
    
if __name__=="__main__":
    name="../Results/Force0000002/8/phase250000.vts"
    #name="../Sailfish/Asymmetric/4/asym200000_0.vti"    
    #contour_movie_vti(name)
    contour_movie(name)
    
    

    