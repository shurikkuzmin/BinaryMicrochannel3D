#include "dynamics_BGK.h"
#include "lattice.h"
#include "solver.h"

void DynamicsBGK::init(int iX,int iY, int iZ)
{
	Lattice * lattice=this->solver->getLattice();

	int zbegin=lattice->getZbegin();
	int zend=lattice->getZend();

	int NX=this->solver->geom->getNX();
	int NY=this->solver->geom->getNY();

	int counter=iZ*NX*NY+iY*NX+iX;

	double gradx_temp=0.0,grady_temp=0.0,gradz_temp=0.0,laplace_temp=0.0;
	for (int k=0;k<NPOP;k++)
	{
		double phase_temp;
		int iZ2=iZ+cz[k];
		int iY2=(iY+cy[k]+NY)%NY;
		int iX2=(iX+cx[k]+NX)%NX;
		if (iZ2==-1)
			phase_temp=lattice->phase_bottom[iY2*NX+iX2];
		else
			if (iZ2==zend-zbegin+1)
				phase_temp=lattice->phase_top[iY2*NX+iX2];
			else
				phase_temp=lattice->phase[iZ2*NX*NY+iY2*NX+iX2];

		gradx_temp+=gradx_stencil[k]*phase_temp;
		grady_temp+=grady_stencil[k]*phase_temp;
		gradz_temp+=gradz_stencil[k]*phase_temp;
		laplace_temp+=laplace_stencil[k]*phase_temp;
	}

	double phase_temp=lattice->phase[counter];
	double dense_temp=lattice->rho[counter];
	double ux_temp=lattice->ux[counter];
	double uy_temp=lattice->uy[counter];
	double uz_temp=lattice->uz[counter];

	double feq;
	double geq;
	double sum=0.0;
	double sum_phase=0.0;
	double phase_square=phase_temp*phase_temp;
	double pressure_bulk=dense_temp/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_temp*laplace_temp;
	double chemical_pot=gammaconst*(aconst*(-phase_temp+phase_temp*phase_temp*phase_temp)-kconst*laplace_temp);

	for (int k=1; k<NPOP; k++)
	{
		feq=weights[k]*(pressure_bulk+dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)
						+1.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+(cz[k]*cz[k]-1.0/3.0)*uz_temp*uz_temp
										 +2.0*(ux_temp*uy_temp*cx[k]*cy[k]+ux_temp*uz_temp*cx[k]*cz[k]+uy_temp*uz_temp*cy[k]*cz[k])))
		+kconst*(wxx[k]*gradx_temp*gradx_temp+wyy[k]*grady_temp*grady_temp+wzz[k]*gradz_temp*gradz_temp
				 +wxy[k]*gradx_temp*grady_temp+wyz[k]*grady_temp*gradz_temp+wzx[k]*gradx_temp*gradz_temp);
		geq=weights[k]*(chemical_pot+phase_temp*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)
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

void DynamicsBGK::updateFields(int iX,int iY,int iZ)
{

   	int NX=this->solver->geom->getNX();
	int NY=this->solver->geom->getNY();
    int counter=iZ*NX*NY+iY*NX+iX;

    Lattice * lattice=this->solver->getLattice();


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

void DynamicsBGK::collide_stream(int iX,int iY,int iZ)
{

    //Calculation of gradients and laplace
    Lattice * lattice=this->solver->getLattice();

    int zbegin=lattice->getZbegin();
    int zend=lattice->getZend();
    int NX=this->solver->geom->getNX();
    int NY=this->solver->geom->getNY();

    int counter=iZ*NX*NY+iY*NX+iX;
	double gradx_temp=0.0,grady_temp=0.0,gradz_temp=0.0,laplace_temp=0.0;

	for (int k=0;k<NPOP;k++)
	{
		int iZ2=iZ+cz[k];
		int iY2=(iY+cy[k]+NY)%NY;
		int iX2=(iX+cx[k]+NX)%NX;
		double phase_temp;
		if (iZ2==-1)
			phase_temp=lattice->phase_bottom[iY2*NX+iX2];
		else
			if (iZ2==zend-zbegin+1)
				phase_temp=lattice->phase_top[iY2*NX+iX2];
			else
				phase_temp=lattice->phase[iZ2*NX*NY+iY2*NX+iX2];

		gradx_temp+=gradx_stencil[k]*phase_temp;
		grady_temp+=grady_stencil[k]*phase_temp;
		gradz_temp+=gradz_stencil[k]*phase_temp;
		laplace_temp+=laplace_stencil[k]*phase_temp;
	}

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


	double phase_square=phase_temp*phase_temp;
	double pressure_bulk=dense_temp/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_temp*laplace_temp;
	double chemical_pot=gammaconst*(aconst*(-phase_temp+phase_temp*phase_temp*phase_temp)-kconst*laplace_temp);

	double sum=0.0;
	double sum_phase=0.0;

	double feq[NPOP], geq[NPOP];
	for (int k=1; k<NPOP; k++)
	{
		feq[k]=weights[k]*(pressure_bulk+dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)
						+1.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+(cz[k]*cz[k]-1.0/3.0)*uz_temp*uz_temp
										 +2.0*(ux_temp*uy_temp*cx[k]*cy[k]+ux_temp*uz_temp*cx[k]*cz[k]+uy_temp*uz_temp*cy[k]*cz[k])))
		+kconst*(wxx[k]*gradx_temp*gradx_temp+wyy[k]*grady_temp*grady_temp+wzz[k]*gradz_temp*gradz_temp
				 +wxy[k]*gradx_temp*grady_temp+wyz[k]*grady_temp*gradz_temp+wzx[k]*gradx_temp*gradz_temp);
		geq[k]=weights[k]*(chemical_pot+phase_temp*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)
						+1.5*phase_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+(cz[k]*cz[k]-1.0/3.0)*uz_temp*uz_temp
										 +2.0*(ux_temp*uy_temp*cx[k]*cy[k]+ux_temp*uz_temp*cx[k]*cz[k]+uy_temp*uz_temp*cy[k]*cz[k])));

		sum+=feq[k];
		sum_phase+=geq[k];
	}

	feq[0]=dense_temp-sum;
	geq[0]=phase_temp-sum_phase;

	//Omega calculation
	double tau_temp=tau_gas+(phase_temp+1.0)/2.0*(tau_liq-tau_gas);
	double omega_rho=1.0/tau_temp;
	double omega_phi=1.0/tau_phi;
	
    //Force calculation
    double sum_force=0.0;
    double force_pop[NPOP];
    for (int k=1;k<NPOP;k++)
    {
        force_pop[k]=(1.0-0.5*omega_rho)*weights[k]*(force_x*((cx[k]-ux_temp)+3.0*cx[k]*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp))+
                force_y*((cy[k]-uy_temp)+3.0*cy[k]*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp))+
                force_z*((cz[k]-uz_temp)+3.0*cz[k]*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)));

        sum_force+=force_pop[k];
    }
    force_pop[0]=-sum_force;

	//Collision procedure
	for (int k=0; k<NPOP; k++)
	{
		//cout<<"Increment="<<-omega*(lattice->f[counter*NPOP+k]-feq[k])<<"\n";
		lattice->f[counter*NPOP+k]+=-omega_rho*(lattice->f[counter*NPOP+k]-feq[k])+force_pop[k];
		lattice->g[counter*NPOP+k]+=-omega_phi*(lattice->g[counter*NPOP+k]-geq[k]);
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

DynamicsBGK::DynamicsBGK(Solver* _solver,ParamsList _params):
	Dynamics(_solver,_params),
	aconst(params("aconst").value<double>()),
	kconst(params("kconst").value<double>()),
	gammaconst(params("gammaconst").value<double>()),
	tau_phi(params("tau_phi").value<double>()),
	tau_liq(params("tau_liq").value<double>()),
	tau_gas(params("tau_gas").value<double>()),
	force_x(params("force_x").value<double>()),
	force_y(params("force_y").value<double>()),
	force_z(params("force_z").value<double>()),
	NPOP(19),
	omega(params("omega").value<double>())
{


	char cx_temp[19]={0,1,-1,0,0,0,0,1,-1,1,-1,0,0,0,0,1,-1,1,-1};
	char cy_temp[19]={0,0,0,1,-1,0,0,1,1,-1,-1,1,-1,1,-1,0,0,0,0};
	char cz_temp[19]={0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1};



	double wxx_temp[19]={0.0,5.0/12.0,5.0/12.0,-1.0/3.0,-1.0/3.0,-1.0/3.0,-1.0/3.0,
		-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,
		-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0};
	double wyy_temp[19]={0.0,-1.0/3.0,-1.0/3.0,5.0/12.0,5.0/12.0,-1.0/3.0,-1.0/3.0,
		-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,
		1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0};
	double wzz_temp[19]={0.0,-1.0/3.0,-1.0/3.0,-1.0/3.0,-1.0/3.0,5.0/12.0,5.0/12.0,
		1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,
		-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0};
	double wxy_temp[19]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		1.0/4.0,-1.0/4.0,-1.0/4.0,1.0/4.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	double wyz_temp[19]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,
		1.0/4.0,-1.0/4.0,-1.0/4.0,1.0/4.0,
		0.0,0.0,0.0,0.0};
	double wzx_temp[19]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		1.0/4.0,-1.0/4.0,-1.0/4.0,1.0/4.0};
	double weights_temp[19]={0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,
		1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,
		1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0};

	for (int counter=0;counter<19;counter++)
	{
		cx[counter]=cx_temp[counter];
		cy[counter]=cy_temp[counter];
		cz[counter]=cz_temp[counter];

		wxx[counter]=wxx_temp[counter];
		wyy[counter]=wyy_temp[counter];
		wzz[counter]=wzz_temp[counter];
		wxy[counter]=wxy_temp[counter];
		wyz[counter]=wyz_temp[counter];
		wzx[counter]=wzx_temp[counter];
		weights[counter]=weights_temp[counter];
	}

	//Parameters for stencils


    double gradx_stencil_temp[19] = {0.0, 1.0/6.0, -1.0/6.0, 0.0, 0.0, 0.0, 0.0, 1.0/12.0, -1.0/12.0, 1.0/12.0,
                                -1.0/12.0, 0.0, 0.0, 0.0, 0.0, 1.0/12.0, -1.0/12.0, 1.0/12.0, -1.0/12.0};

    double grady_stencil_temp[19] = {0.0, 0.0, 0.0, 1.0/6.0, -1.0/6.0, 0.0, 0.0, 1.0/12.0, 1.0/12.0, -1.0/12.0,
                                -1.0/12.0, 1.0/12.0, -1.0/12.0, 1.0/12.0, -1.0/12.0, 0.0, 0.0, 0.0, 0.0};

    double gradz_stencil_temp[19] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0/6.0, -1.0/6.0, 0.0, 0.0, 0.0, 0.0, 1.0/12.0, 1.0/12.0,
                                -1.0/12.0, -1.0/12.0, 1.0/12.0, 1.0/12.0, -1.0/12.0, -1.0/12.0};

    double laplace_stencil_temp[19]={-4.0, 1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/6.0, 1.0/6.0, 1.0/6.0,
                                1.0/6.0, 1.0/6.0,1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0};


	for(int counter=0;counter<19;counter++)
	{
		laplace_stencil[counter]=laplace_stencil_temp[counter];
		gradx_stencil[counter]=gradx_stencil_temp[counter];
		grady_stencil[counter]=grady_stencil_temp[counter];
		gradz_stencil[counter]=gradz_stencil_temp[counter];
	}
}

