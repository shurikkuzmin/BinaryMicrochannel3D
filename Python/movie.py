#! /usr/bin/python
import vtk
import numpy
def contour_movie(name):
    gridreader = vtk.vtkXMLStructuredGridReader()
    gridreader.SetFileName(name)
    gridreader.Update()
    
    grid  = gridreader.GetOutput()
    data  = grid.GetPointData()
    points=grid.GetPoints()
    dims  =grid.GetDimensions()
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
    im_filter=vtk.vtkWindowToImageFilter()
    im_filter.SetInput(renWin)
    #jw->Delete();
    #filter->Delete();
    writer=vtk.vtkJPEGWriter()
    writer.SetFileName("blah.jpg")
    writer.SetInput(im_filter.GetOutput())    
    #writer.SetInputConnection(outline.GetOutputPort())
    writer.Write()

    iren.Start()
       
    
   
    
    
    
if __name__=="__main__":
    name="../Results/Force0000002/8/phase250000.vts"
    contour_movie(name)

    
    

    