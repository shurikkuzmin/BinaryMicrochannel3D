#!/usr/bin/python
import os
import subprocess
import pylab
import numpy
import math
from enthought.tvtk.api import tvtk
from enthought.mayavi import mlab

 
def Run_Simulations():
    print os.getcwd()
    
    force_init=1e-6;
    #os.mkdir("temp")
    #for i in range(1,11):
    #    dir_temp=str(i)
    #    os.mkdir(dir_temp)
    #    subprocess.call(['cp','asymmetric.py',dir_temp+"/"])
    #    os.chdir(dir_temp)
    #    force=force_init*i
    #    subprocess.call(['./asymmetric.py','--force='+str(force),'--batch','--every=50000','--max_iters=200001','--output=asym','--output_format=vtk'])
    #    os.chdir("..")
    for i in range(15,105,5):
        dir_temp=str(i)
        os.mkdir(dir_temp)
        subprocess.call(['cp','asymmetric.py',dir_temp+"/"])
        os.chdir(dir_temp)
        force=force_init*i
        subprocess.call(['./asymmetric.py','--force='+str(force),'--batch','--every=50000','--max_iters=200001','--output=asym','--output_format=vtk'])
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
    #Show_Phase_Slice(phase_numpy,175)
    dims=phase_numpy.shape



    center=phase_numpy[dims[0]/2,dims[1]/2,:]
    z1 = numpy.min(numpy.where(center < 0.0))
    z2 = numpy.max(numpy.where(center < 0.0))
    if z1==0:
        z2=numpy.min(numpy.where(center>0.0))
        z1=numpy.max(numpy.where(center>0.0))+dims[2]
    print z1,z2
    
    mid =((z1+z2)/2)%dims[2]
    print mid
    slice=phase_numpy[:,:,mid]
    shape=slice.shape
    prof_axis=slice[shape[0]/2,:]
    prof_diag=numpy.diag(slice)
    
    axis_zero=(1.0-Get_Zero(prof_axis)/(0.5*(shape[1]-2)))
    diag_zero=(math.sqrt(2)-Get_Zero(prof_diag)/(0.5*(shape[1]-2)))
    
    print "Velx",velx_numpy[dims[0]/2,dims[1]/2,mid]
    print "Vely",vely_numpy[dims[0]/2,dims[1]/2,mid]
    print "Velz",velz_numpy[dims[0]/2,dims[1]/2,mid]
    
    #Calculation of the capillary number
    capillary=velx_numpy[dims[0]/2,dims[1]/2,mid]*2.0/3.0/sqrt(8.0*k*a/9.0)    
    pylab.figure()
    pylab.imshow(phase_numpy[:,:,mid])
    pylab.figure()
    pylab.plot(prof_axis)
    pylab.figure()
    pylab.plot(prof_diag)
    #mlab.show()
    
def Analyze_Consequence():
    print os.getcwd()
    for i in range(1,3):
        dir_temp=str(i)
        os.chdir(dir_temp)
        Analyze_Phase("asym200000_0.vti")
        
        
        #Calculate the capillary number
        os.chdir("..")
    pylab.show()    

def Analyze_Simulations():
    from numpy import genfromtxt
    print os.getcwd()
    capillary_theor=[0.03,0.05,0.08,0.1,0.2,0.4,0.6,0.8]
    capillary_str=["3","5","8","10","20","40","60","80"]
    width_theor=[0.04,0.06,0.08,0.1,0.12,0.13,0.15,0.16]

    exam=[2000,2200,2600,2500,1100,950,600,350]
    good=[0,1,2,3,4,5,6,7]
    #style=["bo","rH","c<","y>","bs","g^"]
    #color=["b","r","c","y","b","g"]
    #style_diff=["b-","r:","c-.","y--","b^","g<"]
    
    #fig=pylab.figure()
    widths=[]
    velocities=[]
    
    #read the giavedoni data (not precise though)
    giavedoni=genfromtxt("planarcasesolution.csv",delimiter=',',dtype=float)[1:]


    for i in range(0, len(capillary_theor)):
        dir_temp="Results/"+capillary_str[i]
        #ratio=float(ny[i]-2)/100.0
        os.chdir(dir_temp)
        #os.chdir("Force")
        #pylab.figure()
        name="capillary200000.npz"
        array=numpy.load(name)
        prof=array['phi'][:,exam[i]]
        #x=numpy.arange(0.0,float(ny[i]))/ratio
        #pylab.imshow(array['phi'])
    
        if i in good:
            #pylab.plot(prof)
            widths.append(Get_Zero(prof))
            vel=array['v'][0]
            vel_prof=vel[:,exam[i]]
            velocities.append(vel_prof[len(vel_prof)/2])
        
        #pylab.plot(x[0:20],prof[0:20],color[i],linewidth=3)
        #pylab.plot(x,prof,style_diff[i],markersize=10,linewidth=3)
        #pylab.figure()
        #pylab.plot(array[:, 520*i])
        #pylab.savefig("grid_phase_prof_"+str(49*i)+".eps", dpi=300)
        #Get_Zero(prof)
        #extrapolator=UnivariateSpline(array[0:(49*i+2)/2, 600*i], numpy.arange(0, (49*i+2)/2),  k=2)
        #print extrapolator(0)

        os.chdir("../..")
    
    fig=pylab.figure()
    capillaries=numpy.array(velocities)*(2.0/3.0)/math.sqrt(8.0*0.04*0.04/9.0)
    print "Widths=",widths
    print "Capillaries=",capillaries
    print "Velocities",velocities
    
    pylab.loglog(giavedoni[:,0],giavedoni[:,1]/2.0,"bD-",linewidth=3,markersize=10)
    pylab.loglog(capillary_theor,width_theor,"ys--",linewidth=3,markersize=10)
    pylab.loglog(capillaries,widths,"go-",linewidth=3,markersize=10)

    pylab.xlim(0.02,1.5)
    pylab.ylim(ymin=0.01)
    numpy.savetxt("capillary.dat",zip(capillaries,widths))
    
    fig.subplots_adjust(left=0.15,bottom=0.15)  
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.xlabel(r'''$Ca$''',fontsize=30)
    pylab.ylabel(r'''$\delta$''',fontsize=30)
    
    #labels=[r'''$H_{eff}='''+str(value-2)+r'''$''' for value in ny]
    pylab.legend(["Giavedoni","Heil","Simulations"],loc=4)
    #pylab.xlim(xmax=15)
    pylab.savefig("capillaries_comparison.eps",format="EPS",dpi=300)

def Analyze_Bubble():
    from numpy import genfromtxt
    print os.getcwd()
    capillary_theor=numpy.array([0.03,0.05,0.08,0.1,0.2,0.4,0.6,0.8])
    capillary_str=numpy.array(["3","5","8","10","20","40","60","80"])
    width_theor=[0.04,0.06,0.08,0.1,0.12,0.13,0.15,0.16]

    exam=[2000,2200,2600,2500,1100,950,600,350]
    good=[0,1,2,3,4,5,6,7]
    good_value=[1,2,4]
    style=["b-.","r--","c."]
    labels=[r'''$Ca='''+ca+r'''$''' for ca in numpy.array(capillary_theor[good_value],dtype=str)]
    #color=["b","r","c","y","b","g"]
    #style_diff=["b-","r:","c-.","y--","b^","g<"]
    
    #fig=pylab.figure()
    low=[1509,1865,381]
    high=[2554,2900,1427]
    widths=[]
    velocities=[]
    fig=pylab.figure()
    for counter,value in enumerate(good_value):
        dir_temp="Results/"+capillary_str[value]
        #ratio=float(ny[i]-2)/100.0
        os.chdir(dir_temp)
        #os.chdir("Force")
        #pylab.figure()
        name="capillary200000.npz"
        array=numpy.load(name)
        thicknesses=[]
        for coor in range(low[counter],high[counter]):
            prof=array['phi'][:,coor]
            thicknesses.append(Get_Zero(prof))
        #prof=array['phi'][:,exam[i]]
        #x=numpy.arange(0.0,float(ny[i]))/ratio
        
        #pylab.imshow(array['phi'])
    
        #if i in good:
            #pylab.plot(prof)
        #    widths.append(Get_Zero(prof))
        #    vel=array['v'][0]
        #    vel_prof=vel[:,exam[i]]
        #    velocities.append(vel_prof[len(vel_prof)/2])
        
        #pylab.plot(x[0:20],prof[0:20],color[i],linewidth=3)
        #pylab.plot(x,prof,style_diff[i],markersize=10,linewidth=3)
        #pylab.figure()
        #pylab.plot(array[:, 520*i])
        #pylab.savefig("grid_phase_prof_"+str(49*i)+".eps", dpi=300)
        #Get_Zero(prof)
        #extrapolator=UnivariateSpline(array[0:(49*i+2)/2, 600*i], numpy.arange(0, (49*i+2)/2),  k=2)
        #print extrapolator(0)
        
        #pylab.figure()
        delta_x=15.0/(len(thicknesses)-1)
        x=delta_x*numpy.arange(0,len(thicknesses))*float(len(thicknesses))/3000.0
        pylab.plot(x,thicknesses,style[counter],linewidth=3)

        os.chdir("../..")
    
    pylab.ylim(ymin=0.0,ymax=0.5)
    pylab.xlim(xmax=5.1)

    leg=pylab.legend(labels)
    legtext = leg.get_texts() # all the text.Text instance in the legend
    for text in legtext:
        text.set_fontsize(30) 
    #set(ltext, fontsize='large') # the legend text fontsize

    fig.subplots_adjust(left=0.15,bottom=0.15)  
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.xlabel(r'''$x$''',fontsize=30)
    pylab.ylabel(r'''$\delta$''',fontsize=30)
    
    #labels=[r'''$H_{eff}='''+str(value-2)+r'''$''' for value in ny]
    #pylab.legend(["Simulations","Giavedoni","Heil"],loc=4)
    #pylab.xlim(xmax=15)
    pylab.savefig("bubble_length.eps",format="EPS",dpi=300)



def Analyze_Velocities():
    print os.getcwd()
    ny=[102,127,152,177,202,227]
    nx=[1501,1876,2251,2626,3001,3376]

    for i in range(0, len(ny)):
        dir_temp="Results/"+str(ny[i])
        ratio=float(ny[i]-2)/100.0
        os.chdir(dir_temp)
        #os.chdir("Force")
        name="grid"+str(200000+50000*i)+".npz"
        array=numpy.load(name)
        ux=array['v'][0]
        #print ux.shape
        prof=ux[:,int(1400*ratio)]
        #pylab.figure()
        #pylab.plot(prof)
        print prof[len(prof)/2]
        os.chdir("../..")


if __name__=="__main__":
    
    Analyze_Consequence()
    #Analyze_Phase("100/asym200000_0.vti")
    #Run_Simulations()
    #Analyze_Simulations()    
    #Analyze_Velocities()
    #Analyze_Bubble()
    #pylab.show()
