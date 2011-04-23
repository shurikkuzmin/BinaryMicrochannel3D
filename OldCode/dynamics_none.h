#ifndef DYNAMICS_NONE
#define DYNAMICS_NONE
#include "lattice.h"
#include "solver.h"

class DynamicsNONE:public Dynamics
{
    public:

    explicit DynamicsNONE(Solver * _solver,ParamsList _params):
        Dynamics(_solver,_params),
        phase_wall(params("phase_wall").value<double>()),
        rho_wall(params("rho_wall").value<double>()),
        NPOP(19)
    {
        double weights_temp[19]={0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,
		1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,
		1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0};

        for (int counter=0;counter<19;counter++)
            weights[counter]=weights_temp[counter];
    }

    virtual ~DynamicsNONE()
    {}

    void collide_stream(int iX,int iY,int iZ)
    {
       	Lattice * lattice=this->solver->getLattice();
        int NX=this->solver->geom->getNX();
        int NY=this->solver->geom->getNY();

        int counter=iZ*NX*NY+iY*NX+iX;
        for (int k=0; k<NPOP; k++)
		{
            lattice->f2[counter*NPOP+k]=lattice->f[counter*NPOP+k];
            lattice->g2[counter*NPOP+k]=lattice->g[counter*NPOP+k];
        }
    }

    virtual void updateFields(int iX,int iY,int iZ)
    {
        int NX=this->solver->geom->getNX();
        int NY=this->solver->geom->getNY();
        Lattice * lattice=this->solver->getLattice();

        int counter=iZ*NX*NY+iY*NX+iX;

        double phase_temp=0.0, dense_temp=0.0, ux_temp=0.0, uy_temp=0.0, uz_temp=0.0;
        for (int k=0;k<NPOP;k++)
        {
            dense_temp+=lattice->f[counter*NPOP+k];
            phase_temp+=lattice->g[counter*NPOP+k];
            ux_temp+=lattice->f[counter*NPOP+k]*cx[k];
            uy_temp+=lattice->f[counter*NPOP+k]*cy[k];
            uz_temp+=lattice->f[counter*NPOP+k]*cz[k];
        }

        ux_temp=ux_temp/dense_temp;
        uy_temp=uy_temp/dense_temp;
        uz_temp=uz_temp/dense_temp;

        lattice->phase[counter]=phase_temp;
        lattice->rho[counter]=dense_temp;
        lattice->ux[counter]=ux_temp;
        lattice->uy[counter]=uy_temp;
        lattice->uz[counter]=uz_temp;
    }

    void init(int iX,int iY, int iZ)
    {
       	Lattice * lattice=this->solver->getLattice();
        int NX=this->solver->geom->getNX();
        int NY=this->solver->geom->getNY();
        int counter=iZ*NX*NY+iY*NX+iX;

        double sum=0.0;
        double sum_phase=0.0;

        double feq;
        double geq;

        for (int k=1; k<NPOP; k++)
        {
            feq=weights[k]*rho_wall;
            geq=weights[k]*phase_wall;

            sum+=feq;
            sum_phase+=geq;

            lattice->f[counter*NPOP+k]=feq;
            lattice->g[counter*NPOP+k]=geq;
        }

        lattice->f[counter*NPOP]=rho_wall-sum;
        lattice->g[counter*NPOP]=phase_wall-sum_phase;

    }

    public:
    const double phase_wall;
    const double rho_wall;
    const int NPOP;

    //TODO: Put velocity set and weights in the Descriptor object
    static const char cx[19];
    static const char cy[19];
    static const char cz[19];


    double weights[19];
};

const char DynamicsNONE::cx[19]={0,1,-1,0, 0,0, 0,1,-1, 1,-1,0, 0, 0, 0,1,-1, 1,-1};
const char DynamicsNONE::cy[19]={0,0, 0,1,-1,0, 0,1, 1,-1,-1,1,-1, 1,-1,0, 0, 0, 0};
const char DynamicsNONE::cz[19]={0,0, 0,0, 0,1,-1,0, 0, 0, 0,1, 1,-1,-1,1, 1,-1,-1};





#endif
