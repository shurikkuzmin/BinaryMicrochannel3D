#!/usr/bin/python
import os
import subprocess
import pylab
import numpy
import math
from enthought.tvtk.api import tvtk
from enthought.mayavi import mlab
import matplotlib.pyplot as plt
plt.rcParams["xtick.major.pad"] = 10
plt.rcParams["ytick.major.pad"] = 10
 
def Run_Simulations():
    print os.getcwd()
    
    force_init=30e-6;
    for i in range(1,5):
        dir_temp=str(i)
        os.mkdir(dir_temp)
        subprocess.call(['cp','steady.py',dir_temp+"/"])
        os.chdir(dir_temp)
        force=force_init*i
        subprocess.call(['./steady.py','--force='+str(force_init),'--batch','--every='+str(50000*i),'--max_iters='+str(200000*i+1),
        '--iwidth=6','--output=steady','--output_format=vtk','--bc_wall_grad_order=1'])
        os.chdir("..")

def Get_Zero(prof):
    zero=0
    #pylab.figure()
    #pylab.plot(prof)
    for counter in range(0, len(prof)/2):
        if prof[counter]>=0 and prof[counter+1]<0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    return zero-0.5

def Read_Phase(name):
    print os.getcwd()
    gridreader = tvtk.XMLImageDataReader()
    gridreader.file_name = name
    gridreader.update()

    grid  = gridreader.output
    data  = grid.point_data
    dims  =grid.dimensions
    dims=dims.tolist()
    dims.reverse()
    phase= numpy.array(data.get_array("phi"))
    velocity=numpy.array(data.get_array("v"))
    velx=velocity[:,0]
    vely=velocity[:,1]
    velz=velocity[:,2]
    phase_numpy=phase.reshape(dims)
    velx_numpy =velx.reshape(dims)
    vely_numpy =vely.reshape(dims)
    velz_numpy =velz.reshape(dims)
    return phase_numpy,velx_numpy,vely_numpy,velz_numpy

def Show_Phase_Slice(phase_numpy,slice):
    fig=pylab.figure()
    pylab.imshow(phase_numpy[:,:,slice])
    
    
def Show_Phase(phase_numpy):
    dims=phase_numpy.shape
    
    fig=mlab.figure()
    src = mlab.pipeline.scalar_field(phase_numpy)
 
    mlab.outline()
    mlab.orientation_axes()
    #v= mlab.pipeline.vector_field(velx_numpy,vely_numpy,velz_numpy)
    #vx=mlab.pipeline.scalar_field(velx_numpy)
    #vy=mlab.pipeline.scalar_field(vely_numpy)
    #vz=mlab.pipeline.scalar_field(velz_numpy)
    extract=mlab.pipeline.extract_grid(src)
    extract.set(x_min=1,x_max=dims[0]-2,y_min=1,y_max=dims[1]-2)
    surf = mlab.pipeline.contour_surface(extract)
 
    #mlab.pipeline.image_plane_widget(vx, plane_orientation='x_axes', slice_index=250)
    #mlab.pipeline.vectors(v, mask_points=20, scale_factor=3.)
    #mlab.pipeline.vector_cut_plane(v, mask_points=2, scale_factor=3)

    
def Analyze_Phase(name):
    
    #parameters of the binary liquid model
    k=0.04
    a=0.04
    
    phase_numpy,velx_numpy,vely_numpy,velz_numpy=Read_Phase(name)
    Show_Phase(phase_numpy)
    dims=phase_numpy.shape


    center=phase_numpy[dims[0]/2,dims[1]/2,:]
    z1 = numpy.min(numpy.where(center < 0.0))
    z2 = numpy.max(numpy.where(center < 0.0))
    if z1==0:
        z2=numpy.min(numpy.where(center>0.0))+dims[2]
        z1=numpy.max(numpy.where(center>0.0))
    pylab.figure()
    pylab.plot(center)
    print z1,z2
    
    mid =((z1+z2)/2)%dims[2]
    print mid
    slice=phase_numpy[:,:,mid]
    shape=slice.shape
    prof_axis=slice[shape[0]/2,:]
    prof_diag=numpy.diag(slice)
    
    axis_zero=(1.0-Get_Zero(prof_axis)/(0.5*(shape[1]-2)))
    diag_zero=(math.sqrt(2.0)-math.sqrt(2.0)*Get_Zero(prof_diag)/(0.5*(shape[1]-2)))
    
    print "Velx",velx_numpy[dims[0]/2,dims[1]/2,mid]
    print "Vely",vely_numpy[dims[0]/2,dims[1]/2,mid]
    print "Velz",velz_numpy[dims[0]/2,dims[1]/2,mid]
    
    #Calculation of the capillary number
    capillary=velx_numpy[dims[0]/2,dims[1]/2,mid]*(2.0/3.0)/math.sqrt(8.0*k*a/9.0)    
    
#    pylab.figure()
#    pylab.imshow(phase_numpy[:,:,mid])
#    pylab.figure()
#    pylab.plot(prof_axis)
#    pylab.figure()
#    pylab.plot(prof_diag)


    print axis_zero,diag_zero
    print "Capillary=",capillary

    return axis_zero,diag_zero,capillary
    #mlab.show()
    
def Analyze_Consequence():
    ax_zeros=[]
    diag_zeros=[]
    capillaries=[]
    for i in range(1,5):
        dir_temp=str(i)
        os.chdir(dir_temp)
        print os.getcwd()
        name="steady"+str(2*i)+"00000_0.vti"
        print name
        axis_zero, diag_zero, capillary=Analyze_Phase(name)
        ax_zeros.append(axis_zero)
        diag_zeros.append(diag_zero)
        capillaries.append(capillary)
        #Calculate the capillary number
        os.chdir("..")
    pylab.plot(capillaries,ax_zeros)
    pylab.plot(capillaries,diag_zeros)
    numpy.savetxt("steady.txt",zip(capillaries,ax_zeros,diag_zeros))
    pylab.show()
    mlab.show()
    
def Analyze_NoForce_Consequence():
    ax_zeros=[]
    diag_zeros=[]
    capillaries=[]
    os.chdir("NoForce")
    for i in range(1,5):
        print os.getcwd()
        name="steady"+str(2*i)+"00000_0.vti"
        print name
        axis_zero, diag_zero, capillary=Analyze_Phase(name)
        ax_zeros.append(axis_zero)
        diag_zeros.append(diag_zero)
        capillaries.append(capillary)
        #Calculate the capillary number
        os.chdir("..")
    pylab.plot(capillaries,ax_zeros)
    pylab.plot(capillaries,diag_zeros)
    #numpy.savetxt("steady.txt",zip(capillaries,ax_zeros,diag_zeros))
    pylab.show()
    mlab.show()


    
        

def Draw_Consequence():
    #from numpy import genfromtxt

    arr=numpy.loadtxt("steady.txt")
    capillaries=arr[:,0]
    axis_zeros=arr[:,1]
    diag_zeros=arr[:,2]
    
    #To define later to compare with real stuff
    #capillary_theor=[0.03,0.05,0.08,0.1,0.2,0.4,0.6,0.8]
    #capillary_str=["3","5","8","10","20","40","60","80"]
    #width_theor=[0.04,0.06,0.08,0.1,0.12,0.13,0.15,0.16]

    
    #read the giavedoni data (not precise though)
    #giavedoni=genfromtxt("planarcasesolution.csv",delimiter=',',dtype=float)[1:]

    fig=pylab.figure()
    
    #pylab.loglog(giavedoni[:,0],giavedoni[:,1]/2.0,"bD-",linewidth=3,markersize=10)
    #pylab.loglog(capillary_theor,width_theor,"ys--",linewidth=3,markersize=10)
    #pylab.loglog(capillaries,widths,"go-",linewidth=3,markersize=10)
    pylab.plot(capillaries,axis_zeros,"bD-",linewidth=3,markersize=10)
    pylab.plot(capillaries,diag_zeros,"ys-",linewidth=3,markersize=10)
    
    #pylab.xlim(0.02,1.5)
    #pylab.ylim(ymin=0.01)
    fig.subplots_adjust(left=0.15,bottom=0.15)  
    #pylab.xticks(fontsize=20)
    #pylab.yticks(fontsize=20)
    pylab.xlabel(r'''$Ca$''',fontsize=30)
    pylab.ylabel(r'''$R_{axis},R_{diag}$''',fontsize=30)
    
    #labels=[r'''$H_{eff}='''+str(value-2)+r'''$''' for value in ny]
    pylab.legend([r'''$R_{axis}$''',r'''$R_{diag}$'''])
    #pylab.xlim(xmax=15)
    pylab.savefig("underresolved_capillaries.eps",format="EPS",dpi=300)

def Analyze_Bubble(name):
    #from numpy import genfromtxt
    print os.getcwd()
 
 
    #parameters of the binary liquid model
    k=0.04
    a=0.04
    
    phase_numpy,velx_numpy,vely_numpy,velz_numpy=Read_Phase(name)
    #Show_Phase_Slice(phase_numpy,175)
    dims=phase_numpy.shape

    center=phase_numpy[dims[0]/2,dims[1]/2,:]
    z1 = numpy.min(numpy.where(center < 0.0))
    z2 = numpy.max(numpy.where(center < 0.0))
    if z1==0:
        z2=numpy.min(numpy.where(center>0.0))+dims[2]
        z1=numpy.max(numpy.where(center>0.0))
    print z1,z2
    mid =((z1+z2)/2)%dims[2]
    print mid
 
    
    #print "Velx",velx_numpy[dims[0]/2,dims[1]/2,mid]
    #print "Vely",vely_numpy[dims[0]/2,dims[1]/2,mid]
    #print "Velz",velz_numpy[dims[0]/2,dims[1]/2,mid]
    
    #Calculation of the capillary number
    capillary=velx_numpy[dims[0]/2,dims[1]/2,mid]*(2.0/3.0)/math.sqrt(8.0*k*a/9.0)    
    print capillary
 
    radiuses_axis=[]
    radiuses_diag=[]
    for z in range(z1+1,z2-1):
        coor=z%dims[2]
        slice=phase_numpy[:,:,coor]
        shape=slice.shape

        prof_axis=slice[shape[0]/2,:]
        prof_diag=numpy.diag(slice)
    
        radiuses_axis.append((1.0-Get_Zero(prof_axis)/(0.5*(shape[1]-2))))
        radiuses_diag.append((math.sqrt(2.0)-math.sqrt(2.0)*Get_Zero(prof_diag)/(0.5*(shape[1]-2))))
   
 
    fig=pylab.figure()
  
    pylab.plot(radiuses_axis,"bD-",linewidth=3,markersize=7)
    pylab.plot(radiuses_diag,"ys-",linewidth=3,markersize=7)
    
    #pylab.ylim(ymin=0.0,ymax=0.5)
    #pylab.xlim(xmax=5.1)
    pylab.legend([r'''$R_{axis}$''',r'''$R_{diag}$'''],loc=6)
 
    #leg=pylab.legend(labels)
   # legtext = leg.get_texts() # all the text.Text instance in the legend
   # for text in legtext:
    #    text.set_fontsize(30) 
    #set(ltext, fontsize='large') # the legend text fontsize

    fig.subplots_adjust(left=0.15,bottom=0.15)  
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.xlabel(r'''$x$''',fontsize=30)
    pylab.ylabel(r'''$R_{axis},R_{diag}$''',fontsize=30)
    
    #labels=[r'''$H_{eff}='''+str(value-2)+r'''$''' for value in ny]
    #pylab.legend(["Simulations","Giavedoni","Heil"],loc=4)
    #pylab.xlim(xmax=15)
    pylab.savefig("bubble_length_ca_one.eps",format="EPS",dpi=300)


if __name__=="__main__":
    
    #Run_Simulations()
    #Analyze_Consequence()
    #Analyze_Simulations()    
    #Analyze_Velocities()
    #Analyze_Bubble()
    #Compare_Plots()
    Analyze_NoForce_Consequence()
    #pylab.show()
