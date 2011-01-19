#include "dynamics_BGK.h"
#include "lattice.h"
#include "solver.h"

template<typename D> void DynamicsBGK<D>::init(int iX,int iY, int iZ)
{
	Lattice<D> * lattice=this->solver->getLattice();

	int zbegin=lattice->getZbegin();
	int zend=lattice->getZend();

	int NX=this->solver->geom->getNX();
	int NY=this->solver->geom->getNY();

	int counter=iZ*NX*NY+iY*NX+iX;

	double gradx_temp=0.0,grady_temp=0.0,gradz_temp=0.0,laplace_temp=0.0;
	for (int k=0;k<D::NPOP;k++)
	{
		double phase_temp;
		int iZ2=iZ+D::cz[k];
		int iY2=(iY+D::cy[k]+NY)%NY;
		int iX2=(iX+D::cx[k]+NX)%NX;
		if (iZ2==-1)
			phase_temp=lattice->phase_bottom[iY2*NX+iX2];
		else
			if (iZ2==zend-zbegin+1)
				phase_temp=lattice->phase_top[iY2*NX+iX2];
			else
				phase_temp=lattice->phase[iZ2*NX*NY+iY2*NX+iX2];

		gradx_temp+=D::gradx_stencil[k]*phase_temp;
		grady_temp+=D::grady_stencil[k]*phase_temp;
		gradz_temp+=D::gradz_stencil[k]*phase_temp;
		laplace_temp+=D::laplace_stencil[k]*phase_temp;
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

	for (int k=1; k<D::NPOP; k++)
	{
		feq=D::weights[k]*(pressure_bulk+dense_temp*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp)
						+1.5*dense_temp*((D::cx[k]*D::cx[k]-1.0/3.0)*ux_temp*ux_temp+(D::cy[k]*D::cy[k]-1.0/3.0)*uy_temp*uy_temp+(D::cz[k]*D::cz[k]-1.0/3.0)*uz_temp*uz_temp
										 +2.0*(ux_temp*uy_temp*D::cx[k]*D::cy[k]+ux_temp*uz_temp*D::cx[k]*D::cz[k]+uy_temp*uz_temp*D::cy[k]*D::cz[k])))
		+kconst*(D::wxx[k]*gradx_temp*gradx_temp+D::wyy[k]*grady_temp*grady_temp+D::wzz[k]*gradz_temp*gradz_temp
				 +D::wxy[k]*gradx_temp*grady_temp+D::wyz[k]*grady_temp*gradz_temp+D::wzx[k]*gradx_temp*gradz_temp);
		geq=D::weights[k]*(chemical_pot+phase_temp*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp)
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

template<typename D>
void DynamicsBGK<D>::updateFields(int iX,int iY,int iZ)
{

   	int NX=this->solver->geom->getNX();
	int NY=this->solver->geom->getNY();
    int counter=iZ*NX*NY+iY*NX+iX;

    Lattice<D> * lattice=this->solver->getLattice();


    double phase_temp=0.0, dense_temp=0.0, ux_temp=0.0, uy_temp=0.0, uz_temp=0.0;
	double * f_temp=&lattice->f[counter*D::NPOP];
	double * g_temp=&lattice->g[counter*D::NPOP];

	dense_temp=f_temp[0]+f_temp[1]+f_temp[2]+f_temp[3]+f_temp[4]
                +f_temp[5]+f_temp[6]+f_temp[7]+f_temp[8]+f_temp[9]
                +f_temp[10]+f_temp[11]+f_temp[12]+f_temp[13]+f_temp[14]+
                +f_temp[15]+f_temp[16]+f_temp[17]+f_temp[18];
	phase_temp=g_temp[0]+g_temp[1]+g_temp[2]+g_temp[3]+g_temp[4]
                +g_temp[5]+g_temp[6]+g_temp[7]+g_temp[8]+g_temp[9]
                +g_temp[10]+g_temp[11]+g_temp[12]+g_temp[13]+g_temp[14]+
                +g_temp[15]+g_temp[16]+g_temp[17]+g_temp[18];

    ux_temp=f_temp[1]-f_temp[2]+f_temp[7]-f_temp[8]+f_temp[9]-f_temp[10]+f_temp[15]-f_temp[16]+f_temp[17]-f_temp[18];
	uy_temp=f_temp[3]-f_temp[4]+f_temp[7]+f_temp[8]-f_temp[9]-f_temp[10]+f_temp[11]-f_temp[12]+f_temp[13]-f_temp[14];
	uz_temp=f_temp[5]-f_temp[6]+f_temp[11]+f_temp[12]-f_temp[13]-f_temp[14]+f_temp[15]+f_temp[16]-f_temp[17]-f_temp[18];
	//for (int k=0;k<D::NPOP;k++)
	//{
		//dense_temp+=lattice->f[counter*D::NPOP+k];
		//phase_temp+=lattice->g[counter*D::NPOP+k];
		//ux_temp+=lattice->f[counter*D::NPOP+k]*D::cx[k];
		//uy_temp+=lattice->f[counter*D::NPOP+k]*D::cy[k];
		//uz_temp+=lattice->f[counter*D::NPOP+k]*D::cz[k];
	//}

	ux_temp=ux_temp/dense_temp;
	uy_temp=uy_temp/dense_temp;
	uz_temp=uz_temp/dense_temp;

	lattice->phase[counter]=phase_temp;
	lattice->rho[counter]=dense_temp;
	lattice->ux[counter]=ux_temp;
	lattice->uy[counter]=uy_temp;
	lattice->uz[counter]=uz_temp;

}

template<typename D>
void DynamicsBGK<D>::collide_stream(int iX,int iY,int iZ)
{

    //Calculation of gradients and laplace
    Lattice<D> * lattice=this->solver->getLattice();

    int zbegin=lattice->getZbegin();
    int zend=lattice->getZend();
    int NX=this->solver->geom->getNX();
    int NY=this->solver->geom->getNY();

    int counter=iZ*NX*NY+iY*NX+iX;
	double gradx_temp=0.0,grady_temp=0.0,gradz_temp=0.0,laplace_temp=0.0;

	for (int k=0;k<D::NPOP;k++)
	{
		int iZ2=iZ+D::cz[k];
		int iY2=(iY+D::cy[k]+NY)%NY;
		int iX2=(iX+D::cx[k]+NX)%NX;
		double phase_temp;
		if (iZ2==-1)
			phase_temp=lattice->phase_bottom[iY2*NX+iX2];
		else
			if (iZ2==zend-zbegin+1)
				phase_temp=lattice->phase_top[iY2*NX+iX2];
			else
				phase_temp=lattice->phase[iZ2*NX*NY+iY2*NX+iX2];

		gradx_temp+=D::gradx_stencil[k]*phase_temp;
		grady_temp+=D::grady_stencil[k]*phase_temp;
		gradz_temp+=D::gradz_stencil[k]*phase_temp;
		laplace_temp+=D::laplace_stencil[k]*phase_temp;
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
    double ux_temp_sq=ux_temp*ux_temp;
    double uy_temp_sq=uy_temp*uy_temp;
    double uz_temp_sq=uz_temp*uz_temp;
    double uxuy_temp=ux_temp*uy_temp;
    double uxuz_temp=ux_temp*uz_temp;
    double uyuz_temp=uy_temp*uz_temp;

    double gradx_temp_sq=gradx_temp*gradx_temp;
    double grady_temp_sq=grady_temp*grady_temp;
    double gradz_temp_sq=gradz_temp*gradz_temp;
    double gradxy_temp=gradx_temp*grady_temp;
    double gradxz_temp=gradx_temp*gradz_temp;
    double gradyz_temp=grady_temp*gradz_temp;

	double feq[D::NPOP], geq[D::NPOP];
	for (int k=1; k<D::NPOP; k++)
	{
		feq[k]=D::weights[k]*(pressure_bulk+dense_temp*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp)
						+1.5*dense_temp*(D::qxx[k]*ux_temp_sq+D::qyy[k]*uy_temp_sq+D::qzz[k]*uz_temp_sq
										 +2.0*(uxuy_temp*D::qxy[k]+uxuz_temp*D::qxz[k]+uyuz_temp*D::qyz[k])))
		+kconst*(D::wxx[k]*gradx_temp_sq+D::wyy[k]*grady_temp_sq+D::wzz[k]*gradz_temp_sq
				 +D::wxy[k]*gradxy_temp+D::wyz[k]*gradyz_temp+D::wzx[k]*gradxz_temp);
		geq[k]=D::weights[k]*(chemical_pot+phase_temp*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp)
						+1.5*phase_temp*(D::qxx[k]*ux_temp_sq+D::qyy[k]*uy_temp_sq+D::qzz[k]*uz_temp_sq
										 +2.0*(uxuy_temp*D::qxy[k]+uxuz_temp*D::qxz[k]+uyuz_temp*D::qyz[k])));

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
    double force_pop[D::NPOP];
    for (int k=1;k<D::NPOP;k++)
    {
        force_pop[k]=(1.0-0.5*omega_rho)*D::weights[k]*(force_x*((D::cx[k]-ux_temp)+3.0*D::cx[k]*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp))+
                force_y*((D::cy[k]-uy_temp)+3.0*D::cy[k]*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp))+
                force_z*((D::cz[k]-uz_temp)+3.0*D::cz[k]*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp)));

        sum_force+=force_pop[k];
    }
    force_pop[0]=-sum_force;

	//Collision procedure
	for (int k=0; k<D::NPOP; k++)
	{
		//cout<<"Increment="<<-omega*(lattice->f[counter*NPOP+k]-feq[k])<<"\n";
		lattice->f[counter*D::NPOP+k]+=-omega_rho*(lattice->f[counter*D::NPOP+k]-feq[k])+force_pop[k];
		lattice->g[counter*D::NPOP+k]+=-omega_phi*(lattice->g[counter*D::NPOP+k]-geq[k]);
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
DynamicsBGK<D>::DynamicsBGK(Solver<D>* _solver,ParamsList _params):
	Dynamics<D>(_solver,_params),
	aconst(_params("aconst").value<double>()),
	kconst(_params("kconst").value<double>()),
	gammaconst(_params("gammaconst").value<double>()),
	tau_phi(_params("tau_phi").value<double>()),
	tau_liq(_params("tau_liq").value<double>()),
	tau_gas(_params("tau_gas").value<double>()),
	force_x(_params("force_x").value<double>()),
	force_y(_params("force_y").value<double>()),
	force_z(_params("force_z").value<double>())
{

}

