#include "dynamics_simple_BGK.h"
#include "lattice.h"
#include "solver.h"

void DynamicsSimpleBGK::init(int iX,int iY, int iZ)
{
	Lattice * lattice=this->solver->getLattice();

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

	for (int k=1; k<NPOP; k++)
	{
		feq=weights[k]*(dense_temp/3.0+dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)
						+1.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+(cz[k]*cz[k]-1.0/3.0)*uz_temp*uz_temp
										 +2.0*(ux_temp*uy_temp*cx[k]*cy[k]+ux_temp*uz_temp*cx[k]*cz[k]+uy_temp*uz_temp*cy[k]*cz[k])));
		geq=weights[k]*(phase_temp/3.0+phase_temp*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)
						+1.5*phase_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+(cz[k]*cz[k]-1.0/3.0)*uz_temp*uz_temp
										 +2.0*(ux_temp*uy_temp*cx[k]*cy[k]+ux_temp*uz_temp*cx[k]*cz[k]+uy_temp*uz_temp*cy[k]*cz[k])));
		sum+=feq;
		sum_phase+=geq;
		lattice->f[counter*NPOP+k]=feq;
		lattice->g[counter*NPOP+k]=geq;
	}
	lattice->f[counter*NPOP]=dense_temp-sum;
	lattice->g[counter*NPOP]=phase_temp-sum_phase;

}

void DynamicsSimpleBGK::updateFields(int iX,int iY,int iZ)
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
void DynamicsSimpleBGK::collide_stream(int iX,int iY,int iZ)
{

            //Calculation of gradients and laplace
           	Lattice * lattice=this->solver->getLattice();

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

			double feq[NPOP], geq[NPOP];
			for (int k=1; k<NPOP; k++)
			{
				feq[k]=weights[k]*(dense_temp/3.0+dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)
								+1.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+(cz[k]*cz[k]-1.0/3.0)*uz_temp*uz_temp
												 +2.0*(ux_temp*uy_temp*cx[k]*cy[k]+ux_temp*uz_temp*cx[k]*cz[k]+uy_temp*uz_temp*cy[k]*cz[k])));
				geq[k]=weights[k]*(phase_temp/3.0+phase_temp*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)
								+1.5*phase_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+(cz[k]*cz[k]-1.0/3.0)*uz_temp*uz_temp
												 +2.0*(ux_temp*uy_temp*cx[k]*cy[k]+ux_temp*uz_temp*cx[k]*cz[k]+uy_temp*uz_temp*cy[k]*cz[k])));

				sum+=feq[k];
				sum_phase+=geq[k];
			}

			feq[0]=dense_temp-sum;
			geq[0]=phase_temp-sum_phase;

            //Force calculation
            double sum_force=0.0;
            double force_pop[NPOP];
            for (int k=1;k<NPOP;k++)
            {
                force_pop[k]=(1.0-0.5*omega)*weights[k]*(force_x*((cx[k]-ux_temp)+3.0*cx[k]*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp))+
 					force_y*((cy[k]-uy_temp)+3.0*cy[k]*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp))+
 					force_z*((cz[k]-uz_temp)+3.0*cz[k]*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)));

                sum_force+=force_pop[k];
            }
            force_pop[0]=-sum_force;

			//Collision procedure
			for (int k=0; k<NPOP; k++)
			{
				//cout<<"Increment="<<-omega*(lattice->f[counter*NPOP+k]-feq[k])<<"\n";
				lattice->f[counter*NPOP+k]+=-omega*(lattice->f[counter*NPOP+k]-feq[k])+force_pop[k];
				lattice->g[counter*NPOP+k]+=-omega*(lattice->g[counter*NPOP+k]-geq[k])+force_pop[k];
			}

			//Streaming procedure (separate loop for normal)
			for (int k=0; k<NPOP; k++)
            {
                int iZ2=(iZ+cz[k]+zend-zbegin+1)%(zend-zbegin+1);
                int iY2=(iY+cy[k]+NY)%NY;
                int iX2=(iX+cx[k]+NX)%NX;
                if (lattice->geometry[iZ2*NX*NY+iY2*NX+iX2]!=SolidNode)
                {
                    lattice->f2[(iZ2*NX*NY+iY2*NX+iX2)*NPOP+k]=lattice->f[counter*NPOP+k];
                    lattice->g2[(iZ2*NX*NY+iY2*NX+iX2)*NPOP+k]=lattice->g[counter*NPOP+k];
                }
            }
}

DynamicsSimpleBGK::DynamicsSimpleBGK(Solver* _solver,ParamsList _params):
	Dynamics(_solver,_params),
	force_x(params("force_x").value<double>()),
	force_y(params("force_y").value<double>()),
	force_z(params("force_z").value<double>()),
	NPOP(19),
	omega(params("omega").value<double>())
{

	char cx_temp[19]={0,1,-1,0,0,0,0,1,-1,1,-1,0,0,0,0,1,-1,1,-1};
	char cy_temp[19]={0,0,0,1,-1,0,0,1,1,-1,-1,1,-1,1,-1,0,0,0,0};
	char cz_temp[19]={0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1};

	double weights_temp[19]={0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,
		1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,
		1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0};

	for (int counter=0;counter<19;counter++)
	{
		cx[counter]=cx_temp[counter];
		cy[counter]=cy_temp[counter];
		cz[counter]=cz_temp[counter];
		weights[counter]=weights_temp[counter];
	}

}

