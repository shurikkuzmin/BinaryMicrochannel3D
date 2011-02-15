#!/usr/bin/python

import numpy
import os
import pylab
from numpy import genfromtxt
 
def draw_capillaries():
    dirs=["Force0000002","Force0000002x82","Force0000005x52"]
    style=["<","^",">"]
    fig=pylab.figure()

    for counter,dir in enumerate(dirs):
        os.chdir(dir)
        print os.getcwd()
        
        mat=numpy.loadtxt("capillaries.txt")
        
        #print mat
        pylab.plot(mat[:,0],0.5*(mat[:,1]+mat[:,2]),style[counter],markersize=10,linewidth=2)
        #pylab.plot(mat[:,0],mat[:,2],"+")
        os.chdir("..")

    
    capillary_theor=[0.905,0.986,1.08,1.19,1.33,1.49]
    radiusses=[0.783,0.778,0.772,0.764,0.753,0.747]
   
   
    #style=["bo","rH","c<","y>","bs","g^"]
    #color=["b","r","c","y","b","g"]
    #style_diff=["b-","r:","c-.","y--","b^","g<"]
    
    #fig=pylab.figure()
    widths=[]
    velocities=[]
    
    #read the giavedoni data (not precise though)
    heil=genfromtxt("heil.csv",delimiter=';',dtype=float)

    #print heil
    #print heil.shape
    
    ax=fig.add_subplot(111)
    #capillaries=numpy.array(velocities)*(2.0/3.0)/math.sqrt(8.0*0.04*0.04/9.0)
    
    pylab.semilogx(heil[:,0],heil[:,1],linewidth=2)
    pylab.semilogx(capillary_theor,radiusses,'o-',markersize=7,linewidth=2)
    #pylab.loglog(giavedoni[:,0],giavedoni[:,1]/2.0,"bD-",linewidth=3,markersize=10)
    #pylab.loglog(capillary_theor,width_theor,"ys--",linewidth=3,markersize=10)
    #pylab.loglog(capillaries,widths,"go-",linewidth=3,markersize=10)

    pylab.xlim(xmin=0.1,xmax=5)
    #pylab.ylim(ymin=0.01)
    #numpy.savetxt("capillary.dat",zip(capillaries,widths))
    
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    
    pylab.ylabel(r'''$R_{diag},R_{axes}$''',fontsize=30)
    pylab.xlabel(r'''$Ca$''',fontsize=30)
    
    
    
    #labels=[r'''$H_{eff}='''+str(value-2)+r'''$''' for value in ny]
    leg=pylab.legend(["CPU results","Refined grid","Large body force","Heil","GPU results"],fancybox=True)
    legtext = leg.get_texts() # all the text.Text instance in the legend
    for text in legtext:
        text.set_fontsize(20) 
    fig.subplots_adjust(left=0.17,bottom=0.17) 

    for line in ax.yaxis.get_majorticklines():
        # line is a Line2D instance
        #line.set_color('green')
        #line.set_markersize(25)
        line.set_markeredgewidth(2)

    for line in ax.xaxis.get_majorticklines():
        # line is a Line2D instance
        #line.set_color('green')
        #line.set_markersize(25)
        line.set_markeredgewidth(2)
    
    for line in ax.yaxis.get_minorticklines():
        # line is a Line2D instance
        #line.set_color('green')
        #line.set_markersize(25)
        line.set_markeredgewidth(2)

    for line in ax.xaxis.get_minorticklines():
        # line is a Line2D instance
        #line.set_color('green')
        #line.set_markersize(25)
        line.set_markeredgewidth(2)


    #for line in ax.get_xticklines() + ax.get_yticklines():
    #    line.set_markersize(10)

 
  

    #pylab.xlim(xmax=15)
    
    pylab.savefig("capillaries_comparison.eps",format="EPS",dpi=300)
    pylab.show()

if __name__=="__main__":
    draw_capillaries()
    pylab.show()