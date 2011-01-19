#ifndef LATTICE_H
#define LATTICE_H

#include "dynamics.h"
#include "geometry.h"
#include <vtkXMLStructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkIdList.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>


template<typename D>
class Lattice
{
    public:
    Lattice(int _zbegin,int _zend,int _NX,int _NY,int _NZ,int _NUMLOCAL):
        zbegin(_zbegin),zend(_zend),
        NX(_NX),NY(_NY),NZ(_NZ),
        NUMLOCAL(_NUMLOCAL)
    {
        //Populations
        f = new double[NUMLOCAL*D::NPOP];
        f2= new double[NUMLOCAL*D::NPOP];
        g = new double[NUMLOCAL*D::NPOP];
        g2= new double[NUMLOCAL*D::NPOP];

        //Macro variables
        phase = new double[NUMLOCAL];
        rho = new double[NUMLOCAL];
        ux  = new double[NUMLOCAL];
        uy  = new double[NUMLOCAL];
        uz  = new double[NUMLOCAL];

        //Helpful materials
        phase_top = new double[NX*NY];
        phase_bottom = new double[NX*NY];
        dynamics_list = new Dynamics<D>*[NUMLOCAL];
        layer_top_f = new double[NX*NY*D::NPOP];
		layer_bottom_f = new double[NX*NY*D::NPOP];
		layer_top_g = new double[NX*NY*D::NPOP];
		layer_bottom_g = new double[NX*NY*D::NPOP];

        geom_top=new NodeType[NX*NY];
        geom_bottom=new NodeType[NX*NY];
        geometry=new NodeType[NUMLOCAL];
    }

    ~Lattice()
    {
        //Populations
        delete[] f;
        delete[] f2;
        delete[] g;
        delete[] g2;

        //Macro variables
        delete[] phase;
        delete[] rho;
        delete[] ux;
        delete[] uy;
        delete[] uz;

        //Helpful materials
        delete[] phase_top;
        delete[] phase_bottom;
        delete[] dynamics_list;
        delete[] layer_top_f;
        delete[] layer_bottom_f;
        delete[] layer_top_g;
        delete[] layer_bottom_g;
        delete[] geom_top;
        delete[] geom_bottom;
        delete[] geometry;
    }

    int getZbegin()
    {
        return zbegin;
    }
    int getZend()
    {
        return zend;
    }
    int getNUMLOCAL()
    {
        return NUMLOCAL;
    }

    void putDensity(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,double _density)
    {
        for(int iZ=iZbegin;iZ<=iZend;iZ++)
            for(int iY=iYbegin;iY<=iYend;iY++)
                for(int iX=iXbegin;iX<=iXend;iX++)
                    rho[iZ*NX*NY+iY*NX+iX]=_density;
    }

    void putPhase(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,double _phase)
    {
        for(int iZ=iZbegin;iZ<=iZend;iZ++)
            for(int iY=iYbegin;iY<=iYend;iY++)
                for(int iX=iXbegin;iX<=iXend;iX++)
                    phase[iZ*NX*NY+iY*NX+iX]=_phase;
    }

    void putVelocity(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,double * _velocity)
    {
        for(int iZ=iZbegin;iZ<=iZend;iZ++)
            for(int iY=iYbegin;iY<=iYend;iY++)
                for(int iX=iXbegin;iX<=iXend;iX++)
                {
                    ux[iZ*NX*NY+iY*NX+iX]=_velocity[0];
                    uy[iZ*NX*NY+iY*NX+iX]=_velocity[1];
                    uz[iZ*NX*NY+iY*NX+iX]=_velocity[2];
                }
    }

    void updateMacro()
    {
        for (int counter=0;counter<NUMLOCAL;counter++)
        {
            int iZ=counter/(NX*NY);
            int iY=(counter%(NX*NY))/NX;
            int iX=(counter%(NX*NY))%NX;
            dynamics_list[counter]->updateFields(iX,iY,iZ);
	    }
    }

    //Output files
    void writePhase(std::string name);
    void writeDensity(std::string name);

    void writeTextPhase(std::string name);
    void writeTextDensity(std::string name);


    //Initialization from the initial field
    void init();

    //Collide and stream
    void collide_stream();

    //Prepare the top and the bottom populations to exchange
    void preparePopulations();
    void finishPropagation();

    //Exchange matrices
    void exchangeMatrices();

    private:
        int zbegin;
        int zend;
        const int NX;
        const int NY;
        const int NZ;
        int NUMLOCAL;

    public:
        double * f;
        double * f2;
        double * g;
        double * g2;
        double * phase;
        double * rho;
        double * ux;
        double * uy;
        double * uz;
        double * phase_top;
        double * phase_bottom;
   		double * layer_top_f;
		double * layer_bottom_f;

		double * layer_top_g;
		double * layer_bottom_g;

        NodeType * geometry;
        NodeType * geom_top;
        NodeType * geom_bottom;

        Dynamics<D> ** dynamics_list;
};
#endif
