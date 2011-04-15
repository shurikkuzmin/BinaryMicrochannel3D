#!/usr/bin/python

import numpy
import os
import pylab
from numpy import genfromtxt
import vtk
import math


def draw_velocities_new():
    files=["vel_moderate.txt","vel_velocity.txt","vel_force0000002.txt"]

    fig=pylab.figure()

    for counter,file in enumerate(files):
        
        mat=numpy.loadtxt(file)
        
        #print mat
        pylab.semilogx(mat[:,0],numpy.divide((mat[:,4]-mat[:,5]),mat[:,4]),"kv",markersize=8)
    
    
    #capillary_theor=[0.905,0.986,1.08,1.19,1.33,1.49]
    #radiusses=[0.783,0.778,0.772,0.764,0.753,0.747]
   
    #style=["bo","rH","c<","y>","bs","g^"]
    #color=["b","r","c","y","b","g"]
    #style_diff=["b-","r:","c-.","y--","b^","g<"]
    
    #fig=pylab.figure()
    widths=[]
    velocities=[]
    
    #read the giavedoni data (not precise though)
    wang=genfromtxt("../Engauge/relative_wang.csv",delimiter=',',dtype=float)

    #print heil
    #print heil.shape
    
    
    #ax=fig.add_subplot(111)
    #capillaries=numpy.array(velocities)*(2.0/3.0)/math.sqrt(8.0*0.04*0.04/9.0)
    
    pylab.semilogx(wang[:,0],wang[:,1],"ko",markersize=8,markerfacecolor="white")
    #pylab.semilogx(capillary_theor,radiusses,'o-',markersize=7,linewidth=2)
    #pylab.loglog(giavedoni[:,0],giavedoni[:,1]/2.0,"bD-",linewidth=3,markersize=10)
    #pylab.loglog(capillary_theor,width_theor,"ys--",linewidth=3,markersize=10)
    #pylab.loglog(capillaries,widths,"go-",linewidth=3,markersize=10)

    #pylab.xlim(xmin=0.1,xmax=5)
    #pylab.ylim(ymin=0.01)
    #numpy.savetxt("capillary.dat",zip(capillaries,widths))
    
    #pylab.xticks(fontsize=20)
    #pylab.yticks(fontsize=20)
    
    #pylab.ylabel(r'''$R_{diag},R_{axes}$''',fontsize=30)
    #pylab.xlabel(r'''$Ca$''',fontsize=30)
    
    
    
    #labels=[r'''$H_{eff}='''+str(value-2)+r'''$''' for value in ny]
    #leg=pylab.legend(["CPU results","Refined grid","Large body force","Heil","GPU results"],fancybox=True)
    #legtext = leg.get_texts() # all the text.Text instance in the legend
    #for text in legtext:
    #    text.set_fontsize(20) 
    #fig.subplots_adjust(left=0.17,bottom=0.17) 

    #for line in ax.yaxis.get_majorticklines():
        # line is a Line2D instance
        #line.set_color('green')
        #line.set_markersize(25)
    #    line.set_markeredgewidth(2)

    #for line in ax.xaxis.get_majorticklines():
        # line is a Line2D instance
        #line.set_color('green')
        #line.set_markersize(25)
    #    line.set_markeredgewidth(2)
    
    #for line in ax.yaxis.get_minorticklines():
        # line is a Line2D instance
        #line.set_color('green')
        #line.set_markersize(25)
    #    line.set_markeredgewidth(2)

    #for line in ax.xaxis.get_minorticklines():
        # line is a Line2D instance
        #line.set_color('green')
        #line.set_markersize(25)
     #   line.set_markeredgewidth(2)


    #for line in ax.get_xticklines() + ax.get_yticklines():
    #    line.set_markersize(10)

 
  

    #pylab.xlim(xmax=15)
    
    #pylab.savefig("capillaries_comparison.eps",format="EPS",dpi=300)
    #pylab.show()


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
    reynolds=[]
    velocities=[]
    superficial=[]
    #file_list=["Moderate/1/phase250000.vts","Moderate/3/phase250000.vts","Moderate/5/phase250000.vts","Moderate/7/phase250000.vts","Moderate/9/phase250000.vts"]
               #"SmallCapGrid/8/phase250000.vts"]
    #file_list=["Velocity/"+str(dir)+"/phase240000.vts" for dir in range(1,11)]    
    file_list=["Force0000002/4/phase240000.vts","Force0000002/6/phase250000.vts","Force0000002/8/phase250000.vts","Force0000002/10/phase250000.vts"]
    #file_list=["Force0000002x82/2/phase180000.vts","Force0000002x82/4/phase180000.vts"] #,"Force0000005x52/8/phase250000.vts","Force0000005x52/10/phase250000.vts"]    
    for file in file_list:
        print file
        axis_zero, diag_zero, capillary,reynolds_value,vel_value,super_value=extract_profiles(file)
        ax_zeros.append(axis_zero)
        diag_zeros.append(diag_zero)
        capillaries.append(capillary)
        reynolds.append(reynolds_value)
        velocities.append(vel_value)
	superficial.append(super_value)
    #pylab.plot(capillaries,ax_zeros)
    #pylab.plot(capillaries,diag_zeros)
    numpy.savetxt("vel_force0000002.txt",zip(capillaries,ax_zeros,diag_zeros,reynolds,velocities,superficial))

    
def draw_velocities():
    names=["force0000002_rel","force0000002x82_rel","force0000005x52_rel"]
    style=["<","^",">"]
    fig1=pylab.figure(1)
    #fig2=pylab.figure(2)
    xcoors=[]
    ycoors=[]
    #numpy.hstack
    for counter,name in enumerate(names):
        if counter>=1:
		    continue
        mat=numpy.loadtxt("Capillaries/cap_"+name+".txt")
        xcoors=numpy.hstack((xcoors,numpy.log10(mat[:,0])))
        #ycoors=numpy.hstack((ycoors,numpy.divide(mat[:,3]-mat[:,5],mat[:,5])))
        ycoors=numpy.divide(mat[:,4]-mat[:,8]/(52.0*52.0),mat[:,4])

        print mat[:,6]+mat[:,7],mat[:,8]
        #print mat.shape
        pylab.figure(1)        
        #pylab.plot(mat[:,0],numpy.divide(mat[:,3]-mat[:,5],mat[:,5]),style[counter],markersize=10,linewidth=2)
        pylab.semilogx(mat[:,0],numpy.divide(mat[:,4]-mat[:,8]/(52.0*52.0),mat[:,4]),style[counter],markersize=10,linewidth=2)
        #pylab.figure(2)
        #pylab.plot(mat[:,0],numpy.divide(mat[:,4]-mat[:,5],mat[:,5]),style[counter],markersize=10,linewidth=2)        
        #pylab.plot(mat[:,0],mat[:,2],"+")

    #Fit a line, ``y = mx + c``, through some noisy data-points:

    #>>> x = np.array([0, 1, 2, 3])
    #>>> y = np.array([-1, 0.2, 0.9, 2.1])

    ##By examining the coefficients, we see that the line should have a
    #gradient of roughly 1 and cuts the y-axis at more-or-less -1.

    #We can rewrite the line equation as ``y = Ap``, where ``A = [[x 1]]``
    #and ``p = [[m], [c]]``.  Now use `lstsq` to solve for `p`:

    A = numpy.vstack([xcoors, numpy.ones(len(xcoors))]).T
    #>>> A
    #array([[ 0.,  1.],
    #       [ 1.,  1.],
    #       [ 2.,  1.],
    #       [ 3.,  1.]])

    a, b = numpy.linalg.lstsq(A, ycoors)[0]
    print a, b
    #1.0 -0.95

    #Plot the data along with the fitted line:

    #>>> import matplotlib.pyplot as plt
    #>>> plt.plot(x, y, 'o', label='Original data', markersize=10)
    pylab.figure(1)
        
    pylab.semilogx(numpy.power(len(xcoors)*[10.0],xcoors), a*xcoors + b, 'r', label='Fitted line')
    #>>> plt.legend()
    #>>> plt.show()
    
        
    #capillary_theor=[0.905,0.986,1.08,1.19,1.33,1.49]
    #radiusses=[0.783,0.778,0.772,0.764,0.753,0.747]
   
   
    #style=["bo","rH","c<","y>","bs","g^"]
    #color=["b","r","c","y","b","g"]
    #style_diff=["b-","r:","c-.","y--","b^","g<"]
    
    #fig=pylab.figure()
    #widths=[]
    #velocities=[]
    
    #read the giavedoni data (not precise though)
    #heil=genfromtxt("Capillaries/heil.csv",delimiter=';',dtype=float)

    #print heil
    #print heil.shape
    
    #ax=fig.add_subplot(111)
    #capillaries=numpy.array(velocities)*(2.0/3.0)/math.sqrt(8.0*0.04*0.04/9.0)
    
    #pylab.semilogx(heil[:,0],heil[:,1],linewidth=2)
    #pylab.semilogx(capillary_theor,radiusses,'o-',markersize=7,linewidth=2)
    #pylab.loglog(giavedoni[:,0],giavedoni[:,1]/2.0,"bD-",linewidth=3,markersize=10)
    #pylab.loglog(capillary_theor,width_theor,"ys--",linewidth=3,markersize=10)
    #pylab.loglog(capillaries,widths,"go-",linewidth=3,markersize=10)

    #pylab.xlim(xmin=0.1,xmax=5)
    #pylab.ylim(ymin=0.01)
    #numpy.savetxt("capillary.dat",zip(capillaries,widths))
    
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    
    pylab.ylabel(r'''$W$''',fontsize=30)
    pylab.xlabel(r'''$Ca$''',fontsize=30)
    
    
    
    #labels=[r'''$H_{eff}='''+str(value-2)+r'''$''' for value in ny]
    #leg=pylab.legend(["CPU results","Refined grid","Large body force","Heil","GPU results"],fancybox=True)
    #legtext = leg.get_texts() # all the text.Text instance in the legend
    #for text in legtext:
    #    text.set_fontsize(20) 
    fig1.subplots_adjust(left=0.17,bottom=0.17) 

    #for line in ax.yaxis.get_majorticklines():
        # line is a Line2D instance
        #line.set_color('green')
        #line.set_markersize(25)
    #    line.set_markeredgewidth(2)

    #for line in ax.xaxis.get_majorticklines():
        # line is a Line2D instance
        #line.set_color('green')
        #line.set_markersize(25)
    #    line.set_markeredgewidth(2)
    
    #for line in ax.yaxis.get_minorticklines():
        # line is a Line2D instance
        #line.set_color('green')
        #line.set_markersize(25)
     #   line.set_markeredgewidth(2)

    #for line in ax.xaxis.get_minorticklines():
        # line is a Line2D instance
        #line.set_color('green')
        #line.set_markersize(25)
    #    line.set_markeredgewidth(2)


    #for line in ax.get_xticklines() + ax.get_yticklines():
    #    line.set_markersize(10)

 
  

    #pylab.xlim(xmax=15)
    
    #pylab.savefig("capillaries_comparison.eps",format="EPS",dpi=300)
    pylab.show()
    
    
    
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

    slug_mid=((z1+z2+dims[2])/2)%dims[2]
    print "Slug_mid=",slug_mid
    phase_mid=numpy.zeros([dims[0],dims[1]])    
    vel_slug=numpy.zeros([dims[0],dims[1]])    
    #for counter in range(0,points.GetNumberOfPoints()):
    #    coorx,coory,coorz=points.GetPoint(counter)
    #    if coorz==mid:
    #       phase_mid[int(coorx),int(coory)]=phase.GetTuple1(counter)
    for coorx in range(0,dims[0]):
        for coory in range(0,dims[1]):
            counter=mid*dims[0]*dims[1]+coory*dims[0]+coorx
            phase_mid[coorx,coory]= phase.GetTuple1(counter)
            counter=slug_mid*dims[0]*dims[1]+coory*dims[0]+coorx
            vel_slug[coorx,coory]= velocity.GetTuple3(counter)[2]
	    
    
	    
	    
    pylab.figure()
    pylab.imshow(phase_mid[1:,1:],cmap="gray",extent=[0.0,1.0,0.0,1.0],origin="lower")
    #Calculation of the capillary number
    capillary=vz[0,z2%dims[2]]*(2.0/3.0)/math.sqrt(8.0*k*a/9.0)    

    number=str(capillary*100)[:2]
    if number[1]==".":
        number=number[:1]
    #pylab.savefig("phase_crossection_ca"+number+".eps")    
    #pylab.imshow(array['phi'],cmap="gray",extent=[0.0,15.0,0.0,1.0])
    #pylab.yticks([0.0,0.5,1.0])

    reynolds=vz[0,z2%dims[2]]*2.0*(dims[0]-2)/(2.0/3.0)
    print "Reynolds=",vz[0,z2%dims[2]]*2.0*(dims[0]-2)/(2.0/3.0)
    prof_axis=phase_mid[0,1:]
    prof_diag=numpy.diag(phase_mid[1:,1:])

    print Get_Zero(prof_axis)
    print Get_Zero(prof_diag)
    axis_zero=Get_Zero(prof_axis)/(dims[0]-2.0)
    diag_zero=math.sqrt(2.0)*Get_Zero(prof_diag)/(dims[0]-2.0)
    
    
    print axis_zero,diag_zero
    print "Capillary=",capillary

    #pylab.figure()
    #pylab.plot(phase_mid[1,1:])    
    #pylab.plot(numpy.diag(phase_mid[1:,1:]))
    
    #area=numpy.ones(phase_mid.shape)    
    #liq=numpy.where(phase_mid>0.0)
    #bubble=numpy.where(phase_mid<0.0)
    #vol_bubble=numpy.sum(velz_mid[bubble])
    #vol_liq=numpy.sum(velz_mid[liq])
    #vol_pos=numpy.sum(velz_bulk)
    #area_bubble=numpy.sum(area[bubble])
    #area_liq=numpy.sum(area[liq])
    
    superficial=float(numpy.sum(vel_slug[1:-1,1:-1]))/float((dims[0]-2.5)*(dims[1]-2.5))    
    print "Interface=",vz[0,z2%dims[2]]
    print "Superficial=",superficial
    return axis_zero,diag_zero,capillary,reynolds,vz[0,z2%dims[2]],superficial


    #pylab.show()

def extract_streamlines():
    ax_zeros=[]
    diag_zeros=[]
    capillaries=[]
 
    dirs=["Force0000002/4","Force0000002/8"]
    for dir in dirs:
        os.chdir(dir)
        print os.getcwd()
        file="phase250000.vts"
        
        if dir=="Force0000002/4":
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
    probe.SpatialMatchOn()    
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
    draw_velocities_new()
    #Analyze_Consequence()
    #extract_streamlines()
    #name="Force0000002/4/phase240000.vts"    
    #read_vtk_2d(name)    
    pylab.show()
