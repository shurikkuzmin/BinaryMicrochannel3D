#!/usr/bin/python

import numpy
import os
import pylab
from numpy import genfromtxt
import vtk
import math
 
def draw_capillaries():
    #dirs=["Force0000002","Force0000002x82","Force0000005x52"]
    #style=["kv","^",">"]
    fig=pylab.figure()

    #for counter,dir in enumerate(dirs):
    #    os.chdir(dir)
    #    print os.getcwd()
    #    
    #    mat=numpy.loadtxt("capillaries.txt")
        
    #    #print mat
    #    pylab.plot(mat[:,0],0.5*(mat[:,1]+mat[:,2]),style[counter],markersize=10,linewidth=2)
        #pylab.plot(mat[:,0],mat[:,2],"+")
    #    os.chdir("..")

    files=["cap_moderate.txt","cap_velocity.txt","cap_force0000002.txt","cap_force0000005x52.txt","cap_force0000002x82.txt"]

    #read the giavedoni data (not precise though)
    heil=genfromtxt("../Engauge/wang.csv",delimiter=',',dtype=float)

    ax=fig.add_subplot(111)
    #capillaries=numpy.array(velocities)*(2.0/3.0)/math.sqrt(8.0*0.04*0.04/9.0)
    
    pylab.semilogx(heil[:,0],heil[:,1],"k>",markersize=7,markerfacecolor="white")
    pylab.semilogx(heil[:,0],heil[:,2],"k<",markersize=7,markerfacecolor="white")
   
    res=numpy.zeros(3)
    for file in files:
        mat=numpy.loadtxt(file)        
        pylab.semilogx(mat[:,0],mat[:,1],"ks",markersize=7)
        pylab.semilogx(mat[:,0],mat[:,2],"ko",markersize=7)
        #res=numpy.vstack((res,mat))
    #capillary_theor=[0.905,0.986,1.08,1.19,1.33,1.49]
    #radiusses=[0.783,0.778,0.772,0.764,0.753,0.747]
   
    #style=["bo","rH","c<","y>","bs","g^"]
    #color=["b","r","c","y","b","g"]
    #style_diff=["b-","r:","c-.","y--","b^","g<"]
    
    #fig=pylab.figure()
    widths=[]
    velocities=[]
    
   
    #print heil
    #print heil.shape
    
    #pylab.semilogx(capillary_theor,radiusses,'o-',markersize=7,linewidth=2)
    #pylab.loglog(giavedoni[:,0],giavedoni[:,1]/2.0,"bD-",linewidth=3,markersize=10)
    #pylab.loglog(capillary_theor,width_theor,"ys--",linewidth=3,markersize=10)
    #pylab.loglog(capillaries,widths,"go-",linewidth=3,markersize=10)

    pylab.xlim(xmin=0.01,xmax=1.2)
    #pylab.ylim(ymin=0.01)
    #numpy.savetxt("capillary.dat",zip(capillaries,widths))
    
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    
    pylab.ylabel(r'''$R_{diag},R_{axis}$''',fontsize=30)
    pylab.xlabel(r'''$Ca$''',fontsize=30)
    
    
    
    #labels=[r'''$H_{eff}='''+str(value-2)+r'''$''' for value in ny]
    leg=pylab.legend([r'''$R_{axis,Wang}$''',r'''$R_{diag,Wang}$''',r'''$R_{axis}$''',r'''$R_{diag}$'''],fancybox=True)
    #legtext = leg.get_texts() # all the text.Text instance in the legend
    #for text in legtext:
    #    text.set_fontsize(20) 
    fig.subplots_adjust(left=0.17,bottom=0.17) 

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
    
    pylab.savefig("capillaries_comparison_wang.eps",format="EPS",dpi=300)
    pylab.show()


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
    file_list=[]
    os.chdir("Force0000002/4")
    for root,dirs,files in os.walk(os.getcwd()):
        for file in files:
            if file[0:5]=="phase":   
                file_list.append(file)
    for file in sorted(file_list):
        print os.getcwd()
        print file
        axis_zero, diag_zero, capillary=extract_profiles(file)
        ax_zeros.append(axis_zero)
        diag_zeros.append(diag_zero)
        capillaries.append(capillary)
        #Calculate the capillary number
        #os.chdir("..")
    #pylab.plot(capillaries,ax_zeros)
    #pylab.plot(capillaries,diag_zeros)
    #numpy.savetxt("steady.txt",zip(capillaries,ax_zeros,diag_zeros))
    pylab.show()
    #mlab.show()

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
    numpy.savetxt("phase.txt",phase_numpy)

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
    capillary=vz[0,mid]*(2.0/3.0)/math.sqrt(8.0*k*a/9.0)    
    
    print axis_zero,diag_zero
    print "Capillary=",capillary

    pylab.figure()
    pylab.plot(numpy.diag(phase_mid[1:,1:]))
    
    return axis_zero,diag_zero,capillary

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


if __name__=="__main__":
    draw_capillaries()
    #Analyze_Consequence()
    #extract_streamlines()
    pylab.show()
