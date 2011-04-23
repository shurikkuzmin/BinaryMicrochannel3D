#include "dynamics_simple_special.h"
#include "solver.h"
DynamicsSimpleSpecial::DynamicsSimpleSpecial(Solver * _solver, ParamsList _params_list):
    DynamicsSimpleBGK(_solver,_params_list)
{
	int NX=this->solver->geom->getNX();
	int NY=this->solver->geom->getNY();
    int NZ=this->solver->geom->getNZ();


    iX=this->solver->iterX;
    iY=this->solver->iterY;
    iZ=this->solver->iterZ;

    for(int k=0; k<NPOP; k++)
    {
        int iX2=(iX+Solver::cx[k]+NX)%NX;
        int iY2=(iY+Solver::cy[k]+NY)%NY;
        int iZ2=(iZ+Solver::cz[k]+NZ)%NZ;
        Geometry * geom=this->solver->geom;
        if (geom->getType(iX2,iY2,iZ2)==SolidNode)
        {
            directions.push_back(k);
        }
    }


    int compliment_temp[19]={0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15};
    for (int k=0;k<NPOP;k++)
        compliment[k]=compliment_temp[k];
}

void DynamicsSimpleSpecial::init(int iX,int iY,int iZ)
{
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

void DynamicsSimpleSpecial::collide_stream(int iX,int iY, int iZ)
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

	double feq[NPOP];
	double geq[NPOP];
	double sum=0.0;
	double sum_phase=0.0;

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

	//Streaming procedure
	for (int k=0; k<NPOP; k++)
	{
   		bool flag=false;
   		for (int j=0;j<directions.size();j++)
   		{
            if(k==directions[j])
            {
                flag=true;
				//cout<<"BB\n";
                lattice->f2[counter*NPOP+compliment[k]]=lattice->f[counter*NPOP+k];
                lattice->g2[counter*NPOP+compliment[k]]=lattice->g[counter*NPOP+k];

            }
        }
   		if (!flag)
   		{
            int iZ2=(iZ+cz[k]+zend-zbegin+1)%(zend-zbegin+1);
            int iY2=(iY+cy[k]+NY)%NY;
            int iX2=(iX+cx[k]+NX)%NX;
            lattice->f2[(iZ2*NX*NY+iY2*NX+iX2)*NPOP+k]=lattice->f[counter*NPOP+k];
            lattice->g2[(iZ2*NX*NY+iY2*NX+iX2)*NPOP+k]=lattice->g[counter*NPOP+k];
   		}
    }

}

