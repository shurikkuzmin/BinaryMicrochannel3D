#include "dynamics_special.h"
#include "solver.h"
DynamicsSpecial::DynamicsSpecial(Solver * _solver, ParamsList _params_list):
    DynamicsBGK(_solver,_params_list),
    phase_gradient(this->params("phase_gradient").value<double>())
{
    normx=0.0;
    normy=0.0;
    normz=0.0;

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
            normx+=double(cx[k]);
            normy+=double(cy[k]);
            normz+=double(cz[k]);
            directions.push_back(k);
        }
    }


    double norm=sqrt(normx*normx+normy*normy+normz*normz);
    normx/=norm;
    normy/=norm;
    normz/=norm;

    wall_densities.resize(NPOP);

    int compliment_temp[19]={0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15};
    for (int k=0;k<NPOP;k++)
        compliment[k]=compliment_temp[k];
}

void DynamicsSpecial::update_all_gradients(int iX,int iY, int iZ)
{
    Lattice * lattice=this->solver->getLattice();
	int zbegin=lattice->getZbegin();
	int zend=lattice->getZend();

	int NX=this->solver->geom->getNX();
	int NY=this->solver->geom->getNY();
    int NZ=this->solver->geom->getNZ();

	int counter=iZ*NX*NY+iY*NX+iX;


    for (int k=1;k<7;k++)
    {
        bool flag=false;
        for(int j=0;j<directions.size();j++)
        {
            if (k==directions[j])
            {
                flag=true;
                int sign=1-2*(k%2);
                grad[k]=phase_gradient*sign;
            }
        }
        if(!flag)
        {
            int sign=1-2*(k%2);

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
                grad[k]=-(phase_temp-lattice->phase[iZ*NX*NY+iY*NX+iX])*sign;
        }
    }
}

void DynamicsSpecial::update_wall_phase(int iX,int iY,int iZ)
{

    //TODO: that's not a proper implementation, especially of wall nodes,
    //      but keep it for a moment.
    Lattice * lattice=this->solver->getLattice();
	int zbegin=lattice->getZbegin();
	int zend=lattice->getZend();

	int NX=this->solver->geom->getNX();
	int NY=this->solver->geom->getNY();

	int counter=iZ*NX*NY+iY*NX+iX;


    double tol=0.00001;

    int dirx;
    int diry;
    int dirz;

    double phasex=0.0;
    double phasey=0.0;
    double phasez=0.0;

    for(int k=0;k<NPOP;k++)
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

        wall_densities[k]=phase_temp;
    }

    if (abs(normx)>tol)
    {
        int signx = (normx > 0) - (normx < 0);
        dirx=-signx;
        int iX2=(iX+dirx+NX)%NX;
        phasex=lattice->phase[iZ*NX*NY+iY*NX+iX2];

        for (int k=0;k<directions.size();k++)
        {
            if (cx[directions[k]]==-dirx)
            {
                wall_densities[directions[k]]=phasex;
            }
        }

    }

    if (abs(normy)>tol)
    {
        int signy = (normy > 0) - (normy < 0);
        diry=-signy;
        int iY2=(iY+diry+NY)%NY;
        phasey=lattice->phase[iZ*NX*NY+iY2*NX+iX];

        for (int k=0;k<directions.size();k++)
        {
            if (cy[directions[k]]==-diry)
            {
                wall_densities[directions[k]]=phasey;
            }
        }

    }

    if (abs(normz)>tol)
    {
        int signz = (normz > 0) - (normz < 0);
        dirz=-signz;
        int iZ2=iZ+dirz;

        if (iZ2==-1)
        {
			phasez=lattice->phase_bottom[iY*NX+iX];
        }
		else
			if (iZ2==zend-zbegin+1)
			{
				phasez=lattice->phase_top[iY*NX+iX];
 			}
			else
			{
				phasez=lattice->phase[iZ2*NX*NY+iY*NX+iX];
			}
        for (int k=0;k<directions.size();k++)
        {
            if (cz[directions[k]]==-dirz)
            {
                //cout<<"The density in the wall with the "<<k<<"directions is "<<wall_densities[directions[k]]<<"\n";
				//cout<<"The gradient is here\n";
                wall_densities[directions[k]]=phasez-2.0*phase_gradient;
            }
        }
    }


}

void DynamicsSpecial::init(int iX,int iY,int iZ)
{
	Lattice * lattice=this->solver->getLattice();

	int zbegin=lattice->getZbegin();
	int zend=lattice->getZend();

	int NX=this->solver->geom->getNX();
	int NY=this->solver->geom->getNY();

	int counter=iZ*NX*NY+iY*NX+iX;

    //Old modification
    //update_wall_phase(iX,iY,iZ);


    //Update gradients
    update_all_gradients(iX,iY,iZ);

	double gradx_temp=0.0,grady_temp=0.0,gradz_temp=0.0,laplace_temp=0.0;

// Old modification calculation through the overall densities
//	for (int k=0;k<19;k++)
//	{
//		gradx_temp+=gradx_stencil[k]*wall_densities[k];
//		grady_temp+=grady_stencil[k]*wall_densities[k];
//		gradz_temp+=gradz_stencil[k]*wall_densities[k];
//		laplace_temp+=laplace_stencil[k]*wall_densities[k];
//	}

    double TOL=0.000001;

    if (abs(normx)>TOL)
        if (normx>0.0)
            gradx_temp=-phase_gradient;
        else
            gradx_temp=phase_gradient;
    else
        gradx_temp=0.5*(grad[1]+grad[2]);


    if (abs(normy)>TOL)
        if (normy>0.0)
            grady_temp=-phase_gradient;
        else
            grady_temp=phase_gradient;
    else
        grady_temp=0.5*(grad[3]+grad[4]);

    if (abs(normz)>0.0)
        if (normz>0.0)
            gradz_temp=-phase_gradient;
        else
            gradz_temp=phase_gradient;
    else
        gradz_temp=0.5*(grad[5]+grad[6]);

    laplace_temp=grad[1]-grad[2]+grad[3]-grad[4]+grad[5]-grad[6];

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

void DynamicsSpecial::collide_stream(int iX,int iY, int iZ)
{
	Lattice * lattice=this->solver->getLattice();

	int zbegin=lattice->getZbegin();
	int zend=lattice->getZend();

	int NX=this->solver->geom->getNX();
	int NY=this->solver->geom->getNY();

	int counter=iZ*NX*NY+iY*NX+iX;

    //Old version
    //update_wall_phase(iX,iY,iZ);

    //Update gradients
    update_all_gradients(iX,iY,iZ);

	double gradx_temp=0.0,grady_temp=0.0,gradz_temp=0.0,laplace_temp=0.0;
//	for (int k=0;k<19;k++)
//	{
//		gradx_temp+=gradx_stencil[k]*wall_densities[k];
//		grady_temp+=grady_stencil[k]*wall_densities[k];
//		gradz_temp+=gradz_stencil[k]*wall_densities[k];
//		laplace_temp+=laplace_stencil[k]*wall_densities[k];
//	}

    double TOL=0.000001;

    if (abs(normx)>TOL)
        if (normx>0.0)
            gradx_temp=-phase_gradient;
        else
            gradx_temp=phase_gradient;
    else
        gradx_temp=0.5*(grad[1]+grad[2]);


    if (abs(normy)>TOL)
        if (normy>0.0)
            grady_temp=-phase_gradient;
        else
            grady_temp=phase_gradient;
    else
        grady_temp=0.5*(grad[3]+grad[4]);

    if (abs(normz)>0.0)
        if (normz>0.0)
            gradz_temp=-phase_gradient;
        else
            gradz_temp=phase_gradient;
    else
        gradz_temp=0.5*(grad[5]+grad[6]);

    laplace_temp=grad[1]-grad[2]+grad[3]-grad[4]+grad[5]-grad[6];


	double phase_temp=lattice->phase[counter];
	double dense_temp=lattice->rho[counter];
	double ux_temp=lattice->ux[counter];
	double uy_temp=lattice->uy[counter];
	double uz_temp=lattice->uz[counter];

	double feq[NPOP];
	double geq[NPOP];
	double sum=0.0;
	double sum_phase=0.0;
	double phase_square=phase_temp*phase_temp;
	double pressure_bulk=dense_temp/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_temp*laplace_temp;
	double chemical_pot=gammaconst*(aconst*(-phase_temp+phase_temp*phase_temp*phase_temp)-kconst*laplace_temp);

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

	double tau_temp=tau_gas+(phase_temp+1.0)/2.0*(tau_liq-tau_gas);
    double omega_temp=1.0/tau_temp;

    //Force calculation
    double sum_force=0.0;
    double force_pop[NPOP];
    for (int k=1;k<NPOP;k++)
    {
        force_pop[k]=(1.0-0.5*omega_temp)*weights[k]*(force_x*((cx[k]-ux_temp)+3.0*cx[k]*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp))+
 					force_y*((cy[k]-uy_temp)+3.0*cy[k]*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp))+
 					force_z*((cz[k]-uz_temp)+3.0*cz[k]*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)));

        sum_force+=force_pop[k];
    }
    force_pop[0]=-sum_force;

	//Collision procedure
	for (int k=0; k<NPOP; k++)
	{
		//cout<<"Increment="<<-omega*(lattice->f[counter*NPOP+k]-feq[k])<<"\n";
		lattice->f[counter*NPOP+k]+=-omega_temp*(lattice->f[counter*NPOP+k]-feq[k])+force_pop[k];
		lattice->g[counter*NPOP+k]+=-omega_temp*(lattice->g[counter*NPOP+k]-geq[k])+force_pop[k];
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

