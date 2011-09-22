#!/usr/bin/python
import numpy
import pylab

def paint_curves():

    wang=numpy.genfromtxt("wang_axial.csv",delimiter=',',dtype=float)
    print wang
    fig=pylab.figure()
    wang2=numpy.array(wang)
    wang2[:,0]=0.0
    print wang
    ind=numpy.where(wang2>1.0)
    print ind
    #ind[1]=ind[1]+1
    wang[ind]=-1.0
    max_length=max(wang[:,0])
    wang[:,0]=max_length-wang[:,0]
    #pylab.plot(wang[:,0],wang[:,3][::-1],"kv",markersize=8) #,markerfacecolor="white")
    pylab.plot(wang[:,0],wang[:,3],"kv",markersize=8) #,markerfacecolor="white")
    pylab.plot(wang[:,0],wang[:,1],"ko",markersize=8) #,markerfacecolor="white")
    pylab.plot(wang[:,0],wang[:,2],"ks",markersize=8) #,markerfacecolor="white")

    #leg=pylab.legend([r'''$Ca=1.0$''',r'''$Ca=0.3$''',r'''$Ca=0.08$'''],loc=3,fancybox=True)
   
    pylab.ylim(ymin=0.0)
    
    #for counter,file in enumerate(files):
        
    #    mat=numpy.loadtxt(file)
        
        #print mat
    #    pylab.semilogx(mat[:,0],numpy.divide((mat[:,4]-mat[:,5]),mat[:,4]),"kv",markersize=8)
    
    
    #capillary_theor=[0.905,0.986,1.08,1.19,1.33,1.49]
    #radiusses=[0.783,0.778,0.772,0.764,0.753,0.747]
   
    #style=["bo","rH","c<","y>","bs","g^"]
    #color=["b","r","c","y","b","g"]
    #style_diff=["b-","r:","c-.","y--","b^","g<"]
    
    
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    
    pylab.ylabel(r'''$R_{axis}$''',fontsize=30)
    pylab.xlabel(r'''$z$''',fontsize=30)
    
    #legtext = leg.get_texts() # all the text.Text instance in the legend
    #for text in legtext:
    #    text.set_fontsize(20) 
    fig.subplots_adjust(left=0.17,bottom=0.17) 



    wang=numpy.genfromtxt("wang_diag.csv",delimiter=',',dtype=float)
    print wang
    fig=pylab.figure()
    #wang2=numpy.array(wang)
    #wang2[:,0]=0.0
    #print wang
    #ind=numpy.where(wang2>1.0)
    #print ind
    #ind[1]=ind[1]+1
    #wang[ind]=-1.0
    pylab.plot(wang[:,0],wang[:,3],"kv",markersize=8) #,markerfacecolor="white")
    pylab.plot(wang[:,0],wang[:,1],"ko",markersize=8) #,markerfacecolor="white")
    pylab.plot(wang[:,0],wang[:,2],"ks",markersize=8) #,markerfacecolor="white")

    leg=pylab.legend([r'''$Ca=1.0$''',r'''$Ca=0.2$''',r'''$Ca=0.06$'''],loc=3,fancybox=True)
   
    pylab.ylim(ymin=0.0,ymax=1.1)
    
    #for counter,file in enumerate(files):
        
    #    mat=numpy.loadtxt(file)
        
        #print mat
    #    pylab.semilogx(mat[:,0],numpy.divide((mat[:,4]-mat[:,5]),mat[:,4]),"kv",markersize=8)
    
    
    #capillary_theor=[0.905,0.986,1.08,1.19,1.33,1.49]
    #radiusses=[0.783,0.778,0.772,0.764,0.753,0.747]
   
    #style=["bo","rH","c<","y>","bs","g^"]
    #color=["b","r","c","y","b","g"]
    #style_diff=["b-","r:","c-.","y--","b^","g<"]
    
    
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    
    pylab.ylabel(r'''$R_{axis}$''',fontsize=30)
    pylab.xlabel(r'''$z$''',fontsize=30)
    
    #legtext = leg.get_texts() # all the text.Text instance in the legend
    #for text in legtext:
    #    text.set_fontsize(20) 
    fig.subplots_adjust(left=0.17,bottom=0.17) 

 
  

    #pylab.xlim(xmax=15)
    
    #pylab.savefig("relative_velocity.eps",format="EPS",dpi=300)
    pylab.show()


if __name__=="__main__":
    paint_curves()
