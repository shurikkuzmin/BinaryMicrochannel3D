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
    k=0.004
    a=0.004
    
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
    for coorx in range(0,dims[0]):
        for coory in range(0,dims[1]):
            counter=mid*dims[0]*dims[1]+coory*dims[0]+coorx
            phase_mid[coorx,coory]= phase.GetTuple1(counter)
    
    pylab.figure()
    pylab.imshow(phase_mid[1:,1:])
    
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
    capillary=vz[0,z2%dims[2]]*(2.0/3.0)/math.sqrt(8.0*k*a/9.0)    
    
    print axis_zero,diag_zero
    print "Capillary=",capillary

    #pylab.figure()
    #pylab.plot(numpy.diag(phase_mid[1:,1:]))
    
    #return axis_zero,diag_zero,capillary

    pylab.show()

def extract_bubble(name):
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


    numpy.savetxt("vz.txt",vz)
    numpy.savetxt("vy.txt",vy)
    numpy.savetxt("phase.txt",phase_numpy[1:,1:])

    #parameters of the binary liquid model
    k=0.004
    a=0.004
    
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
    phase_cur=numpy.zeros([dims[0],dims[1]])    
    for coorx in range(0,dims[0]):
        for coory in range(0,dims[1]):
            counter=mid*dims[0]*dims[1]+coory*dims[0]+coorx
            phase_mid[coorx,coory]= phase.GetTuple1(counter)
    
    pylab.figure()
    pylab.imshow(phase_mid[1:,1:])
    
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
    capillary=vz[0,z2%dims[2]]*(2.0/3.0)/math.sqrt(8.0*k*a/9.0)    
    
    print axis_zero,diag_zero
    print "Capillary=",capillary

    rads_axis=[]
    rads_diag=[]
    
    for coorz in range(z1+2,z2-2):
        z=coorz%dims[2]
        
        for coorx in range(0,dims[0]):
            for coory in range(0,dims[1]):
                counter=z*dims[0]*dims[1]+coory*dims[0]+coorx
                phase_cur[coorx,coory]= phase.GetTuple1(counter)
        
        prof_axis=phase_cur[0,1:]
        prof_diag=numpy.diag(phase_cur[1:,1:])
  
        axis_zero=Get_Zero(prof_axis)/(dims[0]-2.0)
        diag_zero=math.sqrt(2.0)*Get_Zero(prof_diag)/(dims[0]-2.0)
            ##print axis_zero
            ##print diag_zero
        rads_axis.append(axis_zero)
        rads_diag.append(diag_zero) 
    

    fig=pylab.figure()
    #pylab.plot(numpy.diag(phase_mid[1:,1:]))
    pylab.plot(rads_axis,"g-.",linewidth=3,markersize=4)
    pylab.plot(rads_diag,"b-",linewidth=3,markersize=4)
    
    #pylab.ylim(ymin=0.0,ymax=0.5)
    #pylab.xlim(xmax=5.1)
    labels=[r'''$R_{axis}$''',r'''$R_{diag}$''']
    leg=pylab.legend(labels,loc=8)
    legtext = leg.get_texts() # all the text.Text instance in the legend
    for text in legtext:
        text.set_fontsize(30) 
        #set(text, fontsize='large') # the legend text fontsize

    fig.subplots_adjust(left=0.15,bottom=0.15)  
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)

    
    #pylab.figure()
    #pylab.plot(numpy.diag(phase_mid[1:,1:]))
    
    #return axis_zero,diag_zero,capillary

    pylab.show()
    
    
  
if __name__=="__main__":
    #name="../Results/Force0000002/10/phase250000.vts"
    #name="../Results/Capillary01/phase180000.vts"
    name="../Results/LargeCap_thin/phase250000.vts"
    #extract_profiles(name)
    extract_bubble(name)    
    #extract_consequence()
    