#!/usr/bin/python
import vtk
import numpy
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
    capillary=vz[0,mid]*(2.0/3.0)/math.sqrt(8.0*k*a/9.0)    
    
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
    negative=numpy.logical_and(phase_numpy<0.0,y<20)
    x_short=x[::1,::15]
    y_short=y[::1,::15]
    vz_diff_mask=vz_diff #[numpy.where(phase_numpy>0.0)]
    vz_diff_mask[negative]=None
    velx_short=vz_diff_mask[::1,::15]
    vy_mask=vy #[numpy.where(phase_numpy>0.0)]
    vy_mask[negative]=None
    vely_short=vy_mask[::1,::15]

    #pylab.figure()
    #phase_neg=phase_numpy
    #phase_neg[negative]=None
    #pylab.imshow(phase_neg, origin="lower")
    
    pylab.figure(figsize=[31,2.5])
    pylab.quiver(y_short,x_short,velx_short,vely_short,scale=0.6,width=0.0023,headwidth=2)# ,scale=0.1)
    pylab.contour(phase_numpy,[0.0])

def comparison_vortexes():
    dir_name="../Results/Force0000002/"
    for i in range(4,12,2):
        if i==4:
            name=dir_name+str(i)+"/phase240000.vts"
        else:
            name=dir_name+str(i)+"/phase250000.vts"
        extract_profiles(name)
   

if __name__=="__main__":
    name="../Results/Force0000002/4/phase240000.vts"
    #name="../Results/Force0000002/10/phase250000.vts"
    
    #comparison_vortexes()
    extract_profiles(name)
    #extract_bubble(name)
    pylab.show()
 
