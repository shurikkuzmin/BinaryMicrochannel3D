#!/usr/bin/python

import numpy
import os
import pylab
from numpy import genfromtxt
import vtk
import math



def Get_Zero(prof):
    zero=0
    #pylab.figure()
    #pylab.plot(prof)
    for counter in range(0, len(prof)):
        if prof[counter]<0 and prof[counter+1]>=0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    return zero+0.5

def Analyze_Consequence():
    ax_zeros=[]
    diag_zeros=[]
    capillaries=[]
    file_list=["Moderate/3/phase250000.vts","Moderate/7/phase250000.vts","SmallCapGrid/8/phase250000.vts"]
    for file in file_list:
        print file
        axis_zero, diag_zero, capillary=extract_profiles(file)
        ax_zeros.append(axis_zero)
        diag_zeros.append(diag_zero)
        capillaries.append(capillary)
    #pylab.plot(capillaries,ax_zeros)
    #pylab.plot(capillaries,diag_zeros)
    #numpy.savetxt("steady.txt",zip(capillaries,ax_zeros,diag_zeros))
    
def extract_profiles(name):
    
    gridreader = vtk.vtkXMLStructuredGridReader()
    gridreader.SetFileName(name)
    gridreader.Update()
    
    grid  = gridreader.GetOutput()
    data  = grid.GetPointData()
    points=grid.GetPoints()
    dims  =grid.GetDimensions()
    velocity=data.GetArray("Velocity")
    phase=data.GetArray("Phase")
    vz=numpy.zeros([dims[1],dims[2]])
    vy=numpy.zeros([dims[1],dims[2]])
    vx=numpy.zeros([dims[1],dims[2]])
    phase_numpy=numpy.zeros([dims[1],dims[2]])
    
    print vz.shape
    print vy.shape

    for coory in range(0,dims[1]):
        for coorz in range(0,dims[2]):
            counter=coorz*dims[0]*dims[1]+coory*dims[0]
            velx,vely,velz=velocity.GetTuple3(counter)
            vz[coory,coorz]=velz
            vy[coory,coorz]=vely
            vx[coory,coorz]=velx
            phase_numpy[coory,coorz]=phase.GetTuple1(counter)


    #numpy.savetxt("vz.txt",vz)
    #numpy.savetxt("vy.txt",vy)
    #numpy.savetxt("phase.txt",phase_numpy)

    #parameters of the binary liquid model
    k=0.04
    a=0.04
    
    center=phase_numpy[0,:]
    z1 = numpy.min(numpy.where(center < 0.0))
    z2 = numpy.max(numpy.where(center < 0.0))
    if z1==0:
        z2=numpy.min(numpy.where(center>0.0))+dims[2]
        z1=numpy.max(numpy.where(center>0.0))
    print z1,z2
    
    mid =((z1+z2)/2)%dims[2]
    print mid
    
    phase_mid=numpy.zeros([dims[0],dims[1]])    
    #for counter in range(0,points.GetNumberOfPoints()):
    #    coorx,coory,coorz=points.GetPoint(counter)
    #    if coorz==mid:
    #       phase_mid[int(coorx),int(coory)]=phase.GetTuple1(counter)
    for coorx in range(0,dims[0]):
        for coory in range(0,dims[1]):
            counter=mid*dims[0]*dims[1]+coory*dims[0]+coorx
            phase_mid[coorx,coory]= phase.GetTuple1(counter)
    
    #pylab.figure()
    #pylab.imshow(phase_mid[1:,1:],cmap="gray",extent=[0.0,1.0,0.0,1.0],origin="lower")
    #Calculation of the capillary number
    capillary=vz[0,z2%dims[2]]*(2.0/3.0)/math.sqrt(8.0*k*a/9.0)    

    number=str(capillary*100)[:3]
    if number[2]==".":
        number=number[:3]
    pylab.savefig("phase_crossection_ca"+number+".eps")    
    #pylab.imshow(array['phi'],cmap="gray",extent=[0.0,15.0,0.0,1.0])
    #pylab.yticks([0.0,0.5,1.0])

    print "Reynolds=",vz[0,z2%dims[2]]*2.0*(dims[0]-2)/(2.0/3.0)
    prof_axis=phase_mid[0,1:]
    prof_diag=numpy.diag(phase_mid[1:,1:])

    print Get_Zero(prof_axis)
    print Get_Zero(prof_diag)
    axis_zero=Get_Zero(prof_axis)/(dims[0]-2.0)
    diag_zero=math.sqrt(2.0)*Get_Zero(prof_diag)/(dims[0]-2.0)
    
    
    print axis_zero,diag_zero
    print "Capillary=",capillary

    pylab.figure()
    pylab.plot(phase_mid[1,1:])    
    #pylab.plot(numpy.diag(phase_mid[1:,1:]))
    
    return axis_zero,diag_zero,capillary

    #pylab.show()

def extract_streamlines():
    ax_zeros=[]
    diag_zeros=[]
    capillaries=[]
 
    for counter in range(1,11):
        dir="Velocity/"+str(counter)
        os.chdir(dir)
        print os.getcwd()
        file="phase240000.vts"
        
        axis_zero, diag_zero, capillary=extract_profiles(file)
        ax_zeros.append(axis_zero)
        diag_zeros.append(diag_zero)
        capillaries.append(capillary)
        #Calculate the capillary number
        os.chdir("../..")
    #pylab.plot(capillaries,ax_zeros)
    #pylab.plot(capillaries,diag_zeros)
    #numpy.savetxt("steady.txt",zip(capillaries,ax_zeros,diag_zeros))
    pylab.show()
    #mlab.show()

def pyngl_streamline(name):
    from PyNGL import Ngl
    
    gridreader = vtk.vtkXMLStructuredGridReader()
    gridreader.SetFileName(name)
    gridreader.Update()
    
    grid  = gridreader.GetOutput()
    data  = grid.GetPointData()
    points=grid.GetPoints()
    dims  =grid.GetDimensions()
    velocity=data.GetArray("Velocity")
    phase=data.GetArray("Phase")
    vz=numpy.zeros([dims[1],dims[2]])
    vy=numpy.zeros([dims[1],dims[2]])
    vx=numpy.zeros([dims[1],dims[2]])
    phase_numpy=numpy.zeros([dims[1],dims[2]])
    
    print vz.shape
    print vy.shape

    for coory in range(0,dims[1]):
        for coorz in range(0,dims[2]):
            counter=coorz*dims[0]*dims[1]+coory*dims[0]
            velx,vely,velz=velocity.GetTuple3(counter)
            vz[coory,coorz]=velz
            vy[coory,coorz]=vely
            vx[coory,coorz]=velx
            phase_numpy[coory,coorz]=phase.GetTuple1(counter)

    #parameters of the binary liquid model
    k=0.04
    a=0.04
    
    center=phase_numpy[0,:]
    z1 = numpy.min(numpy.where(center < 0.0))
    z2 = numpy.max(numpy.where(center < 0.0))
    if z1==0:
        z2=numpy.min(numpy.where(center>0.0))+dims[2]
        z1=numpy.max(numpy.where(center>0.0))
    print z1,z2
    
    mid =((z1+z2)/2)%dims[2]
    print mid
    
    
    #pylab.figure()
    #pylab.imshow(phase_mid[1:,1:],cmap="gray",extent=[0.0,1.0,0.0,1.0],origin="lower")
    #Calculation of the capillary number
    capillary=vz[0,z2%dims[2]]*(2.0/3.0)/math.sqrt(8.0*k*a/9.0)
    print "Capillary=",capillary
    
    ux=vz-vz[0,z2%dims[2]]        
    uy=vy
    
    #ux_mask=ux[::5,::100]
    #uy_mask=uy[::5,::100]
    #x_mask=x[::5,::100]
    #y_mask=y[::5,::100]
     
     
    wks_type = "eps"
    wks = Ngl.open_wks(wks_type,"stream_ca"+str(100*capillary)[0:2])
    resources = Ngl.Resources()
  
       
    resources.tiMainFont    = "Times-Roman"
    resources.tiMainOn=True
    resources.tiMainString="Ca="+str(capillary)[0:4]
        
    #resources.tiXAxisString = "streamlines"
    resources.vpHeightF = 0.25 # Define height, width, and location of plot.
    resources.vpWidthF  = 3*0.25
    resources.wkPaperSize="A5"
    resources.nglFrame = False
    resources.vfXArray=numpy.linspace(0.0,30.0,len(ux[1,::10]))
    resources.vfYArray=numpy.linspace(0.0,1.0,len(ux[::2,1]))
     
     
    resources2=Ngl.Resources()
    resources2.tiMainFont    = "Times-Roman"
    ##resources2.tiXAxisString = "streamlines"
    resources2.tiMainOn=True
    resources2.tiMainString="Ca="+str(capillary)[0:4]
     
    resources2.wkPaperSize="A5"
    resources2.vpHeightF = 0.25 # Define height, width, and location of plot.
    resources2.vpWidthF  = 3*0.25
    #resources2.nglFrame = False
     
    resources2.cnLineLabelsOn = False   # Turn off contour line labels.
    resources2.cnLinesOn      = True   # Turn off contour lines.
    resources2.cnFillOn       = False    # Turn on contour fill.
    resources2.cnInfoLabelOn   = False 

     
    resources2.cnLevelSelectionMode = "ExplicitLevels"  # Select contour levels. 
    ##resources2.cnMinLevelValF       = 0.0
    ##resources2.cnMaxLevelValF       = 0.001
    ##resources2.cnLevelSpacingF      = 0.0
    resources2.cnLevelCount=1
    resources2.cnLevels=[0.0]
    ##resources2.cnLineThicknesses=[3]
    resources2.cnMonoLineThickness=True
    resources2.cnLineThicknessF=3.0
     
    resources2.lbLabelBarOn=False
    resources2.lbLabelsOn=False
    resources2.sfXArray=numpy.linspace(0.0,30.0,len(ux[1,:]))
    resources2.sfYArray=numpy.linspace(0.0,1.0,len(ux[:,1]))
     
    #plot = Ngl.streamline(wks,uvar[0,::2,::2],vvar[0,::2,::2],resources) 
    #print vz_diff.shape
    #print vy.shape
    #x,y=numpy.mgrid[0:dims[1],0:dims[2]]
     #vx=numpy.sin(x)*numpy.sin(y)
     #vy=numpy.cos(x)*numpy.cos(y)
    #plot=Ngl.streamline(wks,ux[::5,::50],uy[::5,::50],resources)
    plot=Ngl.streamline(wks,ux[::2,::10],uy[::2,::10],resources)
    #Ngl.contour(wks,phase[::5,::50],resources2)        
    Ngl.contour(wks,phase_numpy[:,:],resources2)        
     #plot=Ngl.streamline(wks,vx,vy,resources)
    Ngl.end()
   
    
def read_vtk_2d(name):

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
    data.SetActiveScalars("Phase"); 
    data.SetActiveVectors("Velocity")
    velocity=data.GetArray("Velocity")
    phase = data.GetArray("Phase")
    
    image=vtk.vtkImageData()
    image.SetSpacing(1.0,1.0,1.0)
    image.SetOrigin(0.0,0.0,0.0)
    image.SetDimensions(dims[0],dims[1],dims[2])
    image.GetPointData().SetScalars(phase)
    image.GetPointData().SetVectors(velocity)
    image.Update()
    print "image=",image
    
    extract=vtk.vtkExtractVOI()
    extract.SetInput(image)
    extract.SetVOI(0,0,0,dims[1]-1,0,dims[2]-1)
    extract.Update()
    
    contour=vtk.vtkContourFilter()
    contour.SetInputConnection(extract.GetOutputPort())
    contour.SetValue(0,0.0)
    contour.Update()

    probe=vtk.vtkProbeFilter()
    probe.SetInputConnection(contour.GetOutputPort())    
    probe.SetSource(image)
    #probe.SpatialMatchOn()    
    probe.Update()

    print "Probe=",probe.GetOutput()

    cont=probe.GetOutput()
    vel=cont.GetPointData().GetArray("Velocity")    
    phi=cont.GetPointData().GetArray("Phase")    
    cont_points=cont.GetPoints()
    x_numpy=numpy.zeros(cont_points.GetNumberOfPoints())
    y_numpy=numpy.zeros(cont_points.GetNumberOfPoints())
    z_numpy=numpy.zeros(cont_points.GetNumberOfPoints())    
    
    velx_numpy=numpy.zeros(cont_points.GetNumberOfPoints())
    vely_numpy=numpy.zeros(cont_points.GetNumberOfPoints())
    velz_numpy=numpy.zeros(cont_points.GetNumberOfPoints())
    
    phi_numpy=numpy.zeros(cont_points.GetNumberOfPoints())

    for counter in range(0,cont.GetPoints().GetNumberOfPoints()):
        x,y,z=cont_points.GetPoint(counter)
        x_numpy[counter]=x
        y_numpy[counter]=y
        z_numpy[counter]=z
        velx_numpy[counter]=vel.GetTuple3(counter)[0]
        vely_numpy[counter]=vel.GetTuple3(counter)[1]
        velz_numpy[counter]=vel.GetTuple3(counter)[2]
        phi_numpy[counter]=phi.GetTuple1(counter)
       
    
    
    #Velocity of the interface
    vz=numpy.zeros([dims[1],dims[2]])
    vy=numpy.zeros([dims[1],dims[2]])
    vx=numpy.zeros([dims[1],dims[2]])
    phase_numpy=numpy.zeros([dims[1],dims[2]])
    
    print vz.shape
    print vy.shape

    for coory in range(0,dims[1]):
        for coorz in range(0,dims[2]):
            counter=coorz*dims[0]*dims[1]+coory*dims[0]
            velx,vely,velz=velocity.GetTuple3(counter)
            vz[coory,coorz]=velz
            vy[coory,coorz]=vely
            vx[coory,coorz]=velx
            phase_numpy[coory,coorz]=phase.GetTuple1(counter)

   
    center=phase_numpy[0,:]
    z1 = numpy.min(numpy.where(center < 0.0))
    z2 = numpy.max(numpy.where(center < 0.0))
    if z1==0:
        z2=numpy.min(numpy.where(center>0.0))+dims[2]
        z1=numpy.max(numpy.where(center>0.0))
    print z1,z2
    
    mid =((z1+z2)/2)%dims[2]
    print vz[0,z2%dims[2]]
    print vz[0,((z1+z2)/2)%dims[2]]

    y_numpy=y_numpy/50.0
    z_numpy=z_numpy/50.0
    fig=pylab.figure(figsize=(10,3))
    pylab.plot(z_numpy,y_numpy,"o",markersize=5,color="black")
    pylab.ylim(ymax=1.0)
    #pylab.xlim(xmin=0.1,xmax=5)
    #pylab.ylim(ymin=0.01)
    #numpy.savetxt("capillary.dat",zip(capillaries,widths))
    
    pylab.xticks(fontsize=16)
    pylab.yticks(fontsize=16)
    
    pylab.ylabel(r'''$y$''',fontsize=30)
    pylab.xlabel(r'''$z$''',fontsize=30)
    fig.subplots_adjust(bottom=0.25) 
    pylab.savefig("velocity_interface_contour.eps",dpi=300)

    
    
    #labels=[r'''$H_{eff}='''+str(value-2)+r'''$''' for value in ny]
    #leg=pylab.legend(["CPU results","Refined grid","Large body force","Heil","GPU results"],fancybox=True)
    #legtext = leg.get_texts() # all the text.Text instance in the legend
    #for text in legtext:
    #    text.set_fontsize(20) 

    
    
    #pylab.figure()
    #pylab.plot(z_numpy,phi_numpy,"g+")
    #pylab.figure()
    #pylab.plot(z_numpy,velx_numpy,"+")    
    #pylab.figure()
    #pylab.plot(z_numpy,vely_numpy,"+")    
    fig=pylab.figure(figsize=(10,3))
    pylab.plot(z_numpy,velz_numpy,"o",markersize=5,color="black")
    

    #pylab.plot(z_numpy,y_numpy,"o",markersize=5,color="black")
    
    #pylab.xlim(xmin=0.1,xmax=5)
    #pylab.ylim(ymin=0.01)
    #numpy.savetxt("capillary.dat",zip(capillaries,widths))
    
    pylab.xticks(fontsize=16)
    pylab.yticks(fontsize=16)
    
    pylab.ylabel(r'''$U$''',fontsize=30)
    pylab.xlabel(r'''$z$''',fontsize=30)
    fig.subplots_adjust(bottom=0.25)
    pylab.savefig("velocity_interface_values.eps",dpi=300)
     
    
    #pylab.figure()
    #pylab.plot(vz[0,:],"+")
    
    #pylab.figure()
    #pylab.imshow(phase_numpy)
    #Visualization part 
    #contourMapper=vtk.vtkPolyDataMapper2D()    
    ##contourMapper = vtk.vtkPolyDataMapper()
    ##contourMapper.SetScalarRange(phase.GetRange())
    #contourMapper.SetInputConnection(contour.GetOutputPort())
    ##contourMapper.SetColorMode
 
    ##contourActor = vtk.vtkActor()
    ##contourActor=vtk.vtkActor2D()    
    #contourActor=vtk.vtkImageActor()    
    ##contourActor.SetMapper(contourMapper)
    #contourActor.SetInputConnection(contour.GetOutputPort())
    ##contourActor.SetColor(0.0,0.1,0.2)
    
    
    ## Create the Renderer, RenderWindow, and RenderWindowInteractor
    #ren = vtk.vtkRenderer()
    #renWin = vtk.vtkRenderWindow()
    #renWin.AddRenderer(ren)
    #iren = vtk.vtkRenderWindowInteractor()
    #iren.SetRenderWindow(renWin)

    ## Add the actors to the render; set the background and size
    ##ren.AddActor(streamActor)
    #ren.AddActor(contourActor)
    ##ren.SetBackground(1.0,1.0,1.0)
    #ren.SetBackground(0.1, 0.2, 0.4)
    
    #renWin.SetSize(500, 1500)

    ##ren.ResetCamera()
    ##cam1 = ren.GetActiveCamera()
    ##cam1.SetPosition(1000,0,500)
    ##cam1.SetFocalPoint(0,0,500)

    #iren.Initialize()
    #renWin.Render()
    #iren.Start()


if __name__=="__main__":
    #draw_capillaries()
    #Analyze_Consequence() 
    #name="Velocity/10/phase240000.vts"    
    name="Force0000002/8/phase250000.vts"    
    pyngl_streamline(name)    
    #extract_streamlines()
    #read_vtk_2d(name)    
    #pylab.show()
