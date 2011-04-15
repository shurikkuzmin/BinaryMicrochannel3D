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
    file_list=["Force0000002/10/phase250000.vts","Velocity/1/phase240000.vts","Moderate/3/phase250000.vts"]
    #file_list=["LargeCap/phase250000.vts"]    
    for number,file in enumerate(file_list):
        print file
        print type(number)
        axis_zero, diag_zero, capillary=extract_profiles(file,number)
        ax_zeros.append(axis_zero)
        diag_zeros.append(diag_zero)
        capillaries.append(capillary)
    fig=pylab.figure(1)
    #pylab.ylim(ymin=0.63,ymax=0.7)
    #pylab.xlim(xmax=20)
    labels=[r'''$Ca='''+str(ca)[0:4]+r'''$''' for ca in capillaries]
    ##pylab.plot(rads_axis,"g-.",linewidth=3,markersize=4)
    ##pylab.plot(rads_diag,"b-",linewidth=3,markersize=4)
    
    ##pylab.ylim(ymin=0.0,ymax=0.5)
    ##pylab.xlim(xmax=5.1)
    ##labels=[r'''$R_{axis}$''',r'''$R_{diag}$''']
    #fig=pylab.figure(1)
    pylab.ylabel(r'''$R_{axis}$''',fontsize=20)
    leg=pylab.legend(labels)
    pylab.xticks(fontsize=16)
    pylab.yticks(fontsize=16)
    pylab.xlabel(r'''$z$''',fontsize=20)
    #pylab.title(r'''$Ca='''+str(capillaries[0])[0:4]+r'''$''',fontsize=30)
    #pylab.ylabel(r'''$R_{axis},R_{diag}$''',fontsize=20)
    fig.subplots_adjust(left=0.15,bottom=0.1)  
    pylab.savefig("bubble_rad_axis.eps",dpi=300)
    
    fig=pylab.figure(2)
    pylab.ylabel(r'''$R_{diag}$''',fontsize=20)
    leg=pylab.legend(labels)
    pylab.xticks(fontsize=16)
    pylab.yticks(fontsize=16)
    pylab.xlabel(r'''$z$''',fontsize=20)
    #pylab.ylabel(r'''$y$''',fontsize=20)
       
    fig.subplots_adjust(left=0.15,bottom=0.1)
    pylab.savefig("bubble_rad_diag.eps",dpi=300)
    
    #legtext = leg.get_texts() # all the text.Text instance in the legend
    #for text in legtext:
    #    text.set_fontsize(30) 
    #set(ltext, fontsize='large') # the legend text fontsize

    #pylab.figure(1)
    #leg=str(capillaries)
    #pylab.legend(leg)
    #pylab.plot(capillaries,ax_zeros)
    #pylab.plot(capillaries,diag_zeros)
    #numpy.savetxt("steady.txt",zip(capillaries,ax_zeros,diag_zeros))
    
def extract_profiles(name,filenumber):
    styles=["kv","ko","k^"]
    print "Styles=",styles
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
    
    rad_axis=[]
    rad_diag=[]
    for counter in range(z1+1,z2-1):
        coorz=counter%dims[2]
        
        prof_axis=[]
        prof_diag=[]
        for coory in range(1,dims[1]):
            ind=coorz*dims[0]*dims[1]+coory*dims[0]+1
            prof_axis.append(phase.GetTuple1(ind))
            ind=coorz*dims[0]*dims[1]+coory*dims[0]+coory
            prof_diag.append(phase.GetTuple1(ind))
        
        rad_axis.append(Get_Zero(prof_axis)/(dims[1]-2.0))
        rad_diag.append(math.sqrt(2.0)*Get_Zero(prof_diag)/(dims[1]-2.0))
   
        
    phase_mid=numpy.zeros([dims[0],dims[1]])    
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
    #pylab.savefig("phase_crossection_ca"+number+".eps")    
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

    pylab.figure(1)
    #pylab.plot(phase_mid[1,1:])    
    x=numpy.arange(0,len(rad_axis))*30.0/dims[2]    
    pylab.plot(x,rad_axis,styles[filenumber],markersize=8) #styles[number])
    #pylab.plot(x,rad_axis,"kv",markersize=8) #styles[number])
    pylab.figure(2)    
    pylab.plot(x,rad_diag,styles[filenumber],markersize=8) #styles[number])
    #pylab.plot(x,rad_diag,"k^",markersize=8) #styles[number])
    #pylab.legend([r'''$R_{axis}$''',r'''$R_{diag}$'''])    
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

if __name__=="__main__":
    #draw_capillaries()
    Analyze_Consequence() 
    #name="Force0000005x52/4/phase250000.vts"    
    #extract_profiles(name)
    #name="Force0000002/8/phase250000.vts"    
    #pyngl_streamline(name)    
    #extract_streamlines()
    #read_vtk_2d(name)    
    pylab.show()
