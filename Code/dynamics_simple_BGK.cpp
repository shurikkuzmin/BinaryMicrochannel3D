#include "dynamics_simple_BGK.h"
#include "lattice.h"
#include "solver.h"

template<typename D> void DynamicsSimpleBGK<D>::init(int iX,int iY, int iZ)
{
	Lattice<D> * lattice=this->solver->getLattice();

	int zbegin=lattice->getZbegin();
	int zend=lattice->getZend();

	int NX=this->solver->geom->getNX();
	int NY=this->solver->geom->getNY();

	int counter=iZ*NX*NY+iY*NX+iX;

	double phase_temp=lattice->phase[counter];
	double dense_temp=lattice->rho[counter];
	double ux_temp=lattice->ux[counter];
	double uy_temp=lattice->uy[counter];
	double uz_temp=lattice->uz[counter];

	double feq;
	double geq;
	double sum=0.0;
	double sum_phase=0.0;

	for (int k=1; k<D::NPOP; k++)
	{
		feq=D::weights[k]*(dense_temp/3.0+dense_temp*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp)
						+1.5*dense_temp*((D::cx[k]*D::cx[k]-1.0/3.0)*ux_temp*ux_temp+(D::cy[k]*D::cy[k]-1.0/3.0)*uy_temp*uy_temp+(D::cz[k]*D::cz[k]-1.0/3.0)*uz_temp*uz_temp
										 +2.0*(ux_temp*uy_temp*D::cx[k]*D::cy[k]+ux_temp*uz_temp*D::cx[k]*D::cz[k]+uy_temp*uz_temp*D::cy[k]*D::cz[k])));
		geq=D::weights[k]*(phase_temp/3.0+phase_temp*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp)
						+1.5*phase_temp*((D::cx[k]*D::cx[k]-1.0/3.0)*ux_temp*ux_temp+(D::cy[k]*D::cy[k]-1.0/3.0)*uy_temp*uy_temp+(D::cz[k]*D::cz[k]-1.0/3.0)*uz_temp*uz_temp
										 +2.0*(ux_temp*uy_temp*D::cx[k]*D::cy[k]+ux_temp*uz_temp*D::cx[k]*D::cz[k]+uy_temp*uz_temp*D::cy[k]*D::cz[k])));
		sum+=feq;
		sum_phase+=geq;
		lattice->f[counter*D::NPOP+k]=feq;
		lattice->g[counter*D::NPOP+k]=geq;
	}
	lattice->f[counter*D::NPOP]=dense_temp-sum;
	lattice->g[counter*D::NPOP]=phase_temp-sum_phase;

}

template<typename D> void DynamicsSimpleBGK<D>::updateFields(int iX,int iY,int iZ)
{
    int NX=this->solver->geom->getNX();
    int NY=this->solver->geom->getNY();
    Lattice<D> * lattice=this->solver->getLattice();

    int counter=iZ*NX*NY+iY*NX+iX;

    double phase_temp=0.0, dense_temp=0.0, ux_temp=0.0, uy_temp=0.0, uz_temp=0.0;
	for (int k=0;k<D::NPOP;k++)
	{
		dense_temp+=lattice->f[counter*D::NPOP+k];
		phase_temp+=lattice->g[counter*D::NPOP+k];
		ux_temp+=lattice->f[counter*D::NPOP+k]*D::cx[k];
		uy_temp+=lattice->f[counter*D::NPOP+k]*D::cy[k];
		uz_temp+=lattice->f[counter*D::NPOP+k]*D::cz[k];
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

template<typename D> void DynamicsSimpleBGK<D>::collide_stream(int iX,int iY,int iZ)
{

            //Calculation of gradients and laplace
           	Lattice<D> * lattice=this->solver->getLattice();

            int zbegin=lattice->getZbegin();
            int zend=lattice->getZend();
            int NX=this->solver->geom->getNX();
            int NY=this->solver->geom->getNY();

            int counter=iZ*NX*NY+iY*NX+iX;

 			//Force addition
 			lattice->ux[counter]+=force_x/(2.0*lattice->rho[counter]);
 			lattice->uy[counter]+=force_y/(2.0*lattice->rho[counter]);
 			lattice->uz[counter]+=force_z/(2.0*lattice->rho[counter]);

 			//Reconstruction of equilibrium populations
			double phase_temp=lattice->phase[counter];
			double dense_temp=lattice->rho[counter];
			double ux_temp=lattice->ux[counter];
			double uy_temp=lattice->uy[counter];
			double uz_temp=lattice->uz[counter];

			double sum=0.0;
			double sum_phase=0.0;

			double feq[D::NPOP], geq[D::NPOP];
			for (int k=1; k<D::NPOP; k++)
			{
				feq[k]=D::weights[k]*(dense_temp/3.0+dense_temp*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp)
								+1.5*dense_temp*((D::cx[k]*D::cx[k]-1.0/3.0)*ux_temp*ux_temp+(D::cy[k]*D::cy[k]-1.0/3.0)*uy_temp*uy_temp+(D::cz[k]*D::cz[k]-1.0/3.0)*uz_temp*uz_temp
												 +2.0*(ux_temp*uy_temp*D::cx[k]*D::cy[k]+ux_temp*uz_temp*D::cx[k]*D::cz[k]+uy_temp*uz_temp*D::cy[k]*D::cz[k])));
				geq[k]=D::weights[k]*(phase_temp/3.0+phase_temp*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp)
								+1.5*phase_temp*((D::cx[k]*D::cx[k]-1.0/3.0)*ux_temp*ux_temp+(D::cy[k]*D::cy[k]-1.0/3.0)*uy_temp*uy_temp+(D::cz[k]*D::cz[k]-1.0/3.0)*uz_temp*uz_temp
												 +2.0*(ux_temp*uy_temp*D::cx[k]*D::cy[k]+ux_temp*uz_temp*D::cx[k]*D::cz[k]+uy_temp*uz_temp*D::cy[k]*D::cz[k])));

				sum+=feq[k];
				sum_phase+=geq[k];
			}

			feq[0]=dense_temp-sum;
			geq[0]=phase_temp-sum_phase;

            //Force calculation
            double sum_force=0.0;
            double force_pop[D::NPOP];
            for (int k=1;k<D::NPOP;k++)
            {
                force_pop[k]=(1.0-0.5*omega)*D::weights[k]*(force_x*((D::cx[k]-ux_temp)+3.0*D::cx[k]*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp))+
 					force_y*((D::cy[k]-uy_temp)+3.0*D::cy[k]*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp))+
 					force_z*((D::cz[k]-uz_temp)+3.0*D::cz[k]*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp)));

                sum_force+=force_pop[k];
            }
            force_pop[0]=-sum_force;

			//Collision procedure
			for (int k=0; k<D::NPOP; k++)
			{
				//cout<<"Increment="<<-omega*(lattice->f[counter*NPOP+k]-feq[k])<<"\n";
				lattice->f[counter*D::NPOP+k]+=-omega*(lattice->f[counter*D::NPOP+k]-feq[k])+force_pop[k];
				lattice->g[counter*D::NPOP+k]+=-omega*(lattice->g[counter*D::NPOP+k]-geq[k])+force_pop[k];
			}

			//Streaming procedure (separate loop for normal)
			for (int k=0; k<D::NPOP; k++)
            {
                int iZ2=(iZ+D::cz[k]+zend-zbegin+1)%(zend-zbegin+1);
                int iY2=(iY+D::cy[k]+NY)%NY;
                int iX2=(iX+D::cx[k]+NX)%NX;
                if (lattice->geometry[iZ2*NX*NY+iY2*NX+iX2]!=SolidNode)
                {
                    lattice->f2[(iZ2*NX*NY+iY2*NX+iX2)*D::NPOP+k]=lattice->f[counter*D::NPOP+k];
                    lattice->g2[(iZ2*NX*NY+iY2*NX+iX2)*D::NPOP+k]=lattice->g[counter*D::NPOP+k];
                }
            }
}

template<typename D>
DynamicsSimpleBGK<D>::DynamicsSimpleBGK(Solver<D>* _solver,ParamsList _params):
	Dynamics<D>(_solver,_params),
	force_x(_params("force_x").value<double>()),
	force_y(_params("force_y").value<double>()),
	force_z(_params("force_z").value<double>()),
	omega(_params("omega").value<double>())
{

}

