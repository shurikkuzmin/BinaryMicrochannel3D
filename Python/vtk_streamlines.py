#!/usr/bin/python
import vtk
import numpy
from PyNGL import Ngl

import matplotlib
#matplotlib.use('GTKAgg')
import pylab
import math
import matplotlib.pyplot as plt
plt.rcParams["xtick.major.pad"] = 10
plt.rcParams["ytick.major.pad"] = 10

#plt.rcParams["text.usetex"]=True

def extract_files(name):
    #very long and expensive procedure
    
    gridreader = vtk.vtkXMLStructuredGridReader()
    gridreader.SetFileName(name)
    gridreader.Update()
    
    grid  = gridreader.GetOutput()
    data  = grid.GetPointData()
    points=grid.GetPoints()
    dims  =grid.GetDimensions()
    velocity=data.GetArray("Velocity")
    phase=data.GetArray("Phase")
    vz=numpy.zeros([dims[0],dims[1],dims[2]])
    vy=numpy.zeros([dims[0],dims[1],dims[2]])
    vx=numpy.zeros([dims[0],dims[1],dims[2]])
    phase_numpy=numpy.zeros([dims[0],dims[1],dims[2]])
    
    print vz.shape
    print vy.shape

    for counter in range(0,points.GetNumberOfPoints()):
        coorx,coory,coorz=points.GetPoint(counter)
        velx,vely,velz=velocity.GetTuple3(counter)
        vz[int(coorx),int(coory),int(coorz)]=velz
        vy[int(coorx),int(coory),int(coorz)]=vely
        vx[int(coorx),int(coory),int(coorz)]=velx
        phase_numpy[int(coorx),int(coory),int(coorz)]=phase.GetTuple1(counter)

    numpy.savetxt("vz.txt",vz[0,:,:])
    numpy.savetxt("vy.txt",vy[0,:,:])
    numpy.savetxt("phase.txt",phase_numpy[0,:,:])

def Get_Zero(prof):
    zero=0
    #pylab.figure()
    #pylab.plot(prof)
    for counter in range(0, len(prof)):
        if prof[counter]<0 and prof[counter+1]>=0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    return zero+0.5


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

    #for counter in range(0,points.GetNumberOfPoints()):
    #    coorx,coory,coorz=points.GetPoint(counter)
    #    if coorx==0:
    #        velx,vely,velz=velocity.GetTuple3(counter)
    #       vz[int(coory),int(coorz)]=velz
    #        vy[int(coory),int(coorz)]=vely
    #        vx[int(coory),int(coorz)]=velx
    #       phase_numpy[int(coory),int(coorz)]=phase.GetTuple1(counter)

    for coory in range(0,dims[1]):
        for coorz in range(0,dims[2]):
            counter=coorz*dims[0]*dims[1]+coory*dims[0]
            velx,vely,velz=velocity.GetTuple3(counter)
            vz[coory,coorz]=velz
            vy[coory,coorz]=vely
            vx[coory,coorz]=velx
            phase_numpy[coory,coorz]=phase.GetTuple1(counter)


    numpy.savetxt("vz.txt",vz)
    numpy.savetxt("vy.txt",vy)
     
    numpy.savetxt("phase.txt",phase_numpy[1:,1:])

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
    #pylab.imshow(phase_mid[1:,1:])
    
    prof_axis=phase_mid[0,1:]
    prof_diag=numpy.diag(phase_mid[1:,1:])

    print Get_Zero(prof_axis)
    print Get_Zero(prof_diag)
    axis_zero=Get_Zero(prof_axis)/(dims[0]-2.0)
    diag_zero=math.sqrt(2.0)*Get_Zero(prof_diag)/(dims[0]-2.0)
    
    print "Velx",vx[0,mid]
    print "Vely",vy[0,mid]
    print "Velz",vz[0,mid]
    
    #Calculation of the capillary number
    #capillary=vz[0,mid]*(2.0/3.0)/math.sqrt(8.0*k*a/9.0)    
    capillary=vz[0,z2%dims[2]]*(2.0/3.0)/math.sqrt(8.0*k*a/9.0)
    print axis_zero,diag_zero
    print "Capillary=",capillary
    vel_bulk=vz[0,((z1+z2+dims[2])/2)%dims[2]]
    print "Velocity_bulk=",vel_bulk
    print "Velocity difference=",vz[0,z2%dims[2]]-vel_bulk
    
    vz_diff=vz-vz[0,z2%dims[2]]
    numpy.savetxt("vz_diff.txt",vz-vz[0,z2%dims[2]])
 
    
    #fig=pylab.figure()
    #pylab.title(r'''$Ca='''+str(capillary)[0:4]+r'''$''',size=30)
    #ax1 = fig.add_subplot(111)
    #ax1.plot(phase_mid[0,:],"b-",linewidth=3)
    #ax1.set_ylabel("Phase",color='b',size=30)
    #ax2 = ax1.twinx()
    #ax2.plot(numpy.abs(vy[:,(z2+20)%dims[2]+1]-vy[:,(z2+20)%dims[2]-1]),"r.",linewidth=3)
    #ax2.set_ylabel(r'''$\partial_x v_y$''',color='r',size=30)    
 
    #fig=pylab.figure()
    #pylab.plot(vy[:,(z2+20)%dims[2]])
    #pylab.plot(vz[0,:]-vz[0,z2%dims[2]])
 
    #fig=pylab.figure()
    #pylab.imshow(phase_numpy,origin="lower")
    x,y=numpy.mgrid[0:dims[0],0:dims[2]]
    positive=numpy.where(phase_numpy>0.0)
    negative=numpy.where(phase_numpy<0.0)
    large=numpy.where(x>25)
    x_short=x[::3,::15]
    y_short=y[::3,::15]
    vz_diff_mask=vz_diff #[numpy.where(phase_numpy>0.0)]
    vz_diff_mask[negative]=None
    vz_diff_mask[large]=None
    
    velx_short=vz_diff_mask[::3,::15]
    vy_mask=vy #[numpy.where(phase_numpy>0.0)]
    vy_mask[negative]=None
    vy_mask[large]=None
    
    vely_short=vy_mask[::3,::15]

    #pylab.figure()
    #phase_neg=phase_numpy
    #phase_neg[negative]=None
    #pylab.imshow(phase_neg, origin="lower")
    
    fig=pylab.figure(figsize=[23,2.5])
    fig=fig.subplots_adjust(bottom=0.15,top=0.8)
    pylab.quiver(y_short,x_short,velx_short,vely_short,scale=0.1,width=0.0023,headwidth=2)# ,scale=0.1)
    pylab.title(r'''$Ca='''+str(capillary)[0:4]+r'''$''',size=30)
    pylab.contour(phase_numpy,[0.0],linewidths=[4])
    if capillary>=1.0:
        name_cap=str(capillary*100)[0:3]
    else:
        name_cap=str(capillary*100)[0:2]
    pylab.savefig("vortex_ca"+name_cap+".eps",format="EPS",dpi=300)
    pylab.show()

def extract_profiles_ngl(name):

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

    #for counter in range(0,points.GetNumberOfPoints()):
    #    coorx,coory,coorz=points.GetPoint(counter)
    #    if coorx==0:
    #        velx,vely,velz=velocity.GetTuple3(counter)
    #       vz[int(coory),int(coorz)]=velz
    #        vy[int(coory),int(coorz)]=vely
    #        vx[int(coory),int(coorz)]=velx
    #       phase_numpy[int(coory),int(coorz)]=phase.GetTuple1(counter)

    for coory in range(0,dims[1]):
        for coorz in range(0,dims[2]):
            counter=coorz*dims[0]*dims[1]+coory*dims[0]
            velx,vely,velz=velocity.GetTuple3(counter)
            vz[coory,coorz]=velz
            vy[coory,coorz]=vely
            vx[coory,coorz]=velx
            phase_numpy[coory,coorz]=phase.GetTuple1(counter)


    numpy.savetxt("vz.txt",vz)
    numpy.savetxt("vy.txt",vy)
     
    numpy.savetxt("phase.txt",phase_numpy[1:,1:])

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
    
    
    prof_axis=phase_mid[0,1:]
    prof_diag=numpy.diag(phase_mid[1:,1:])

    print Get_Zero(prof_axis)
    print Get_Zero(prof_diag)
    axis_zero=Get_Zero(prof_axis)/(dims[0]-2.0)
    diag_zero=math.sqrt(2.0)*Get_Zero(prof_diag)/(dims[0]-2.0)
    
    print "Velx",vx[0,mid]
    print "Vely",vy[0,mid]
    print "Velz",vz[0,mid]
    
    #Calculation of the capillary number
    capillary=vz[0,mid]*(2.0/3.0)/math.sqrt(8.0*k*a/9.0)    
    
    print axis_zero,diag_zero
    print "Capillary=",capillary
    vel_bulk=vz[0,((z1+z2+dims[2])/2)%dims[2]]
    print "Velocity_bulk=",vel_bulk
    print "Velocity difference=",vz[0,z2%dims[2]]-vel_bulk
    
    vz_diff=vz-vz[0,z2%dims[2]]
    numpy.savetxt("vz_diff.txt",vz-vz[0,z2%dims[2]])
 
    
    x,y=numpy.mgrid[0:dims[0],0:dims[2]]
    positive=numpy.where(phase_numpy>0.0)
    negative=numpy.where(phase_numpy<0.0)
    large=numpy.where(x>35)
    x_short=x[::3,::15]
    y_short=y[::3,::15]
    vz_diff_mask=vz_diff #[numpy.where(phase_numpy>0.0)]
    vz_diff_mask[negative]=None
    vz_diff_mask[large]=None
    
    velx_short=vz_diff_mask[::3,::15]
    vy_mask=vy #[numpy.where(phase_numpy>0.0)]
    vy_mask[negative]=None
    vy_mask[large]=None
    
    vely_short=vy_mask[::3,::15]

    wks_type = "ps"
    wks = Ngl.open_wks(wks_type,"vortex_ngl")
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
    #resources.vpHeightF = 0.25 # Define height, width, and location of plot.
    #resources.vpWidthF  = 7.5*0.25
     #plot = Ngl.streamline(wks,uvar[0,::2,::2],vvar[0,::2,::2],resources) 
    print vz_diff.shape
    print vy.shape
    plot=Ngl.streamline(wks,vz_diff[:,::50],vy[:,::50],resources)
    Ngl.end()

def extract_profiles_vtk(name):

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
    
    extract=vtk.vtkExtractVOI()
    extract.SetInput(phase_image)
    extract.SetVOI(0,0,0,dims[1]-1,0,dims[2]-1)
    extract.Update()
     
    #Create a contour    
    #contour=vtk.vtkOutlineFilter()
    contour=vtk.vtkContourFilter()
    contour.SetInputConnection(extract.GetOutputPort())
    #contour.SetInput(phase_image)
    contour.SetValue(0,0.0)
    contour.Update()
   
    #cont_points=contour.GetOutput().GetPoints()
    probe=vtk.vtkProbeFilter()
    probe.SetInputConnection(contour.GetOutputPort())
    probe.SetSource(extract.GetOutput())
    probe.Update()
    
    filt_points=vtk.vtkMaskPoints()
    filt_points.SetInputConnection(probe.GetOutputPort())    
    #filt_points.SetMaximumNumberOfPoints(200)
    #filt_points.SetRandomMode(1)
    filt_points.SetOnRatio(10)
    arrow=vtk.vtkArrowSource()
    glyph=vtk.vtkGlyph3D()
    glyph.SetInputConnection(filt_points.GetOutputPort())
    glyph.SetSourceConnection(arrow.GetOutputPort())
    glyph.SetVectorModeToUseVector()
    glyph.SetScaleModeToScaleByVector()
    glyph.SetScaleFactor(400.0)
    
    glyphMapper=vtk.vtkPolyDataMapper()
    glyphMapper.SetInputConnection(glyph.GetOutputPort())

    glyphActor=vtk.vtkActor()
    glyphActor.SetMapper(glyphMapper)
    glyphActor.GetProperty().SetColor(0.1,0.5,0.2)
    print glyphActor.GetProperty().GetColor()    
    #print "Probe=", probe.GetOutput()
    calc=vtk.vtkArrayCalculator()
    calc.SetInput(probe.GetOutput())
    #calc.AddVectorArrayName("Velocity",1,1,1)    
    calc.AddScalarVariable("Vz", "Velocity", 2)
    calc.SetResultArrayName("Velocity comp 2")
    calc.SetFunction("Vz")    
    #calc.SetResultArrayName("Velocity Magnitude")
    #calc.SetFunction("mag(Velocity)")
    calc.Update()
    print calc.GetOutput()
    
    xyplot = vtk.vtkXYPlotActor()
    vel=probe.GetOutput().GetPointData().GetArray("Velocity")
    
    #exam=vtk.vtkDoubleArray()
    #vel.GetData(0,vel.GetNumberOfTuples()-1,2,2,exam)
    xyplot.AddInput(calc.GetOutput())
    #xyplot.GetPositionCoordinate().SetValue(0.0, 0.67, 0)
    #xyplot.GetPosition2Coordinate().SetValue(1.0, 0.33, 0) #relative to Position
    xyplot.SetXValuesToArcLength()
    #xyplot.SetNumberOfXLabels(6)
    #xyplot.SetTitle("Pressure vs. Arc Length (Zoomed View)") 
    #xyplot.SetXTitle("")
    #xyplot.SetYTitle("P")
    #xyplot.SetXRange(.1, .35)
    #xyplot.SetYRange(.2, .4)
    xyplot.GetProperty().SetColor(0, 0, 0)
    xyplot.GetProperty().SetLineWidth(2)
    # Set text prop color (same color for backward compat with test) 
    # Assign same object to all text props
    tprop = xyplot.GetTitleTextProperty()
    tprop.SetColor(xyplot.GetProperty().GetColor())
    xyplot.SetAxisTitleTextProperty(tprop)
    xyplot.SetAxisLabelTextProperty(tprop)
    
    #glyph=vtk.vtkGlyph3D()
    #glyph.SetInput(contour.GetOutput())
 
    
    contourMapper = vtk.vtkPolyDataMapper()
    #contourMapper.SetScalarRange(phase.GetRange())
    contourMapper.SetInputConnection(contour.GetOutputPort())
    
    sliceMapper = vtk.vtkImageMapper()
    #sliceMapper = vtk.vtkDataSetMapper()
    sliceMapper.SetInput(extract.GetOutput())
    sliceMapper.SetColorLevel(1000)    
    sliceMapper.SetColorWindow(2000)    
    #sliceMapper.SetColorModeToMapScalars()
 
    #print polydata.GetClassName()
    
    contourActor = vtk.vtkActor()
    contourActor.SetMapper(contourMapper)
   
    sliceActor = vtk.vtkActor2D()
    sliceActor.SetMapper(sliceMapper)

    ren = vtk.vtkRenderer()
    #ren.AddActor(contourActor)
    ren.AddActor2D(sliceActor)
    #ren.AddActor(glyphActor)
    ren.SetBackground(0.1, 0.2, 0.4)
    #ren.SetColor(0.1,0.5,0.2)    
    #ren.SetViewport(0, 0, .3, 1)

    ren2=vtk.vtkRenderer()
    ren2.SetBackground(1, 1, 1)
    ren2.SetViewport(0.3, 0.0, 1.0, 1.0)
    ren2.AddActor2D(xyplot)

    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    #renWin.AddRenderer(ren2)
    renWin.SetSize(500, 500)
    
    #win=vtk.vtkWindowToImageFilter()
    #win.SetInput(renWin)
    
    
    
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    
    iren.Initialize()
    renWin.Render()
    iren.Start()
    
    #save=vtk.vtkPostScriptWriter()
    #save.SetInput(win.GetOutput())
    #save.SetFileName("big_window.eps")
    #save.Write()
    
    
    #vz=numpy.zeros([dims[1],dims[2]])
    #vy=numpy.zeros([dims[1],dims[2]])
    #vx=numpy.zeros([dims[1],dims[2]])
    #phase_numpy=numpy.zeros([dims[1],dims[2]])
    
    #print vz.shape
    #print vy.shape

    #for counter in range(0,points.GetNumberOfPoints()):
    #    coorx,coory,coorz=points.GetPoint(counter)
    #    if coorx==0:
    #        velx,vely,velz=velocity.GetTuple3(counter)
    #       vz[int(coory),int(coorz)]=velz
    #        vy[int(coory),int(coorz)]=vely
    #        vx[int(coory),int(coorz)]=velx
    #       phase_numpy[int(coory),int(coorz)]=phase.GetTuple1(counter)

    #for coory in range(0,dims[1]):
    #    for coorz in range(0,dims[2]):
    #        counter=coorz*dims[0]*dims[1]+coory*dims[0]
    #        velx,vely,velz=velocity.GetTuple3(counter)
    #       vz[coory,coorz]=velz
    #        vy[coory,coorz]=vely
    #        vx[coory,coorz]=velx
    #        phase_numpy[coory,coorz]=phase.GetTuple1(counter)


    #parameters of the binary liquid model
    #k=0.04
    #a=0.04
    
    #center=phase_numpy[0,:]
    #z1 = numpy.min(numpy.where(center < 0.0))
    #z2 = numpy.max(numpy.where(center < 0.0))
    #if z1==0:
    #    z2=numpy.min(numpy.where(center>0.0))+dims[2]
    #    z1=numpy.max(numpy.where(center>0.0))
    #print z1,z2
    
    #mid =((z1+z2)/2)%dims[2]
    #print mid
    
    #phase_mid=numpy.zeros([dims[0],dims[1]])    
    
    #for counter in range(0,points.GetNumberOfPoints()):
    #    coorx,coory,coorz=points.GetPoint(counter)
    #    if coorz==mid:
    #       phase_mid[int(coorx),int(coory)]=phase.GetTuple1(counter)
    #for coorx in range(0,dims[0]):
    #    for coory in range(0,dims[1]):
    #        counter=mid*dims[0]*dims[1]+coory*dims[0]+coorx
    #        phase_mid[coorx,coory]= phase.GetTuple1(counter)
    
    #pylab.figure()
    #pylab.imshow(phase_mid[1:,1:])
    
    #prof_axis=phase_mid[0,1:]
    #prof_diag=numpy.diag(phase_mid[1:,1:])

    #print Get_Zero(prof_axis)
    #print Get_Zero(prof_diag)
    #axis_zero=Get_Zero(prof_axis)/(dims[0]-2.0)
    #diag_zero=math.sqrt(2.0)*Get_Zero(prof_diag)/(dims[0]-2.0)
    
    #print "Velx",vx[0,mid]
    #print "Vely",vy[0,mid]
    #print "Velz",vz[0,mid]
    
    #Calculation of the capillary number
    #capillary=vz[0,mid]*(2.0/3.0)/math.sqrt(8.0*k*a/9.0)    
    
    #print axis_zero,diag_zero
    #print "Capillary=",capillary
    #vel_bulk=vz[0,((z1+z2+dims[2])/2)%dims[2]]
    #print "Velocity_bulk=",vel_bulk
    #print "Velocity difference=",vz[0,z2%dims[2]]-vel_bulk
    
    #vz_diff=vz-vz[0,z2%dims[2]]
    #numpy.savetxt("vz_diff.txt",vz-vz[0,z2%dims[2]])
 
     
def comparison_vortexes():
    dir_name="../Results/Force0000002/"
    for i in range(4,12,2):
        if i==4:
            name=dir_name+str(i)+"/phase240000.vts"
        else:
            name=dir_name+str(i)+"/phase250000.vts"
        extract_profiles(name)
   

if __name__=="__main__":
    name="../Results/Force0000002/8/phase250000.vts"
    #name="../Results/Force0000002/10/phase250000.vts"
    
    #comparison_vortexes()
    extract_profiles(name)
    #extract_profiles_vtk(name)
    #extract_bubble(name)
    #pylab.show()
 
