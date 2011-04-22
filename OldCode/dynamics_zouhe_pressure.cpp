#include "dynamics_zouhe_pressure.h"
#include "solver.h"

DynamicsZouHePressure::DynamicsZouHePressure(Solver * _solver,ParamsList _params_list):
    DynamicsBGK(_solver,_params_list),
    normx(params("normx").value<int>()),
    normy(params("normy").value<int>()),
    normz(params("normz").value<int>()),
    rho_press(params("rho_press").value<double>()),
    phase_press(params("phase_press").value<double>())
{

    int compliment_temp[19]={0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15};
    for (int k=0;k<NPOP;k++)
        compliment[k]=compliment_temp[k];

     for(int iCoor=0;iCoor<19;iCoor++)
    {
        //Determination of the projection
        int dot_product=cx[iCoor]*normx+cy[iCoor]*normy+cz[iCoor]*normz;
        if (dot_product>0)
        {
            unknown.push_back(iCoor);
            positive.push_back(iCoor);
            if ( ! ((cx[iCoor]==normx)&&(cy[iCoor]==normy)&&(cz[iCoor]==normz)) )
            {
                diagonal.push_back(iCoor);
            }
        }
        else
        {
            known.push_back(iCoor);
            if (dot_product<0)
                negative.push_back(iCoor);
            else
                neutral.push_back(iCoor);
        }
    }

}

void DynamicsZouHePressure::updateFields(int iX,int iY,int iZ)
{
  	int NX=this->solver->geom->getNX();
	int NY=this->solver->geom->getNY();
    Lattice * lattice=this->solver->getLattice();


    int counter=iZ*NX*NY+iY*NX+iX;
    lattice->rho[counter]=rho_press;
    lattice->phase[counter]=phase_press;
}

void DynamicsZouHePressure::collide_stream(int iX,int iY,int iZ)
{
    int NX=this->solver->geom->getNX();
    int NY=this->solver->geom->getNY();
    Lattice * lattice=this->solver->getLattice();
   	int zbegin=lattice->getZbegin();
	int zend=lattice->getZend();



    int counter=iZ*NX*NY+iY*NX+iX;
    double denseminus=0.0;
    double phaseminus=0.0;
    double denseneutral=0.0;
    double phaseneutral=0.0;

    for(int iCount=0;iCount<negative.size();iCount++)
    {
        denseminus+=lattice->f[counter*NPOP+negative[iCount]];
        phaseminus+=lattice->g[counter*NPOP+negative[iCount]];
    }

    for(int iCount=0;iCount<neutral.size();iCount++)
    {
        denseneutral+=lattice->f[counter*NPOP+neutral[iCount]];
        phaseneutral+=lattice->g[counter*NPOP+neutral[iCount]];
    }

    double velocity_perp=(2.0*denseminus+denseneutral-rho_press)/rho_press;
    double velocity_tang=0.0;

    //Transform velocities to real velocities
    double norma=sqrt(normx*normx+normy*normy+normz*normz);
    double ux_temp=velocity_perp*normx/norma;
    double uy_temp=velocity_perp*normy/norma;
    double uz_temp=velocity_perp*normz/norma;

	//Calculation of gradients and so on
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


	double phase_square=phase_press*phase_press;
	double pressure_bulk=rho_press/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_press*laplace_temp;
	double chemical_pot=gammaconst*(aconst*(-phase_press+phase_press*phase_press*phase_press)-kconst*laplace_temp);

    for(int iCount=0;iCount<unknown.size();iCount++)
    {
        int pos=unknown[iCount];

        double eqpos=weights[pos]*(pressure_bulk+rho_press*(cx[pos]*ux_temp+cy[pos]*uy_temp+cz[pos]*uz_temp)
						+1.5*rho_press*((cx[pos]*cx[pos]-1.0/3.0)*ux_temp*ux_temp
                                         +(cy[pos]*cy[pos]-1.0/3.0)*uy_temp*uy_temp
                                         +(cz[pos]*cz[pos]-1.0/3.0)*uz_temp*uz_temp
										 +2.0*(ux_temp*uy_temp*cx[pos]*cy[pos]
										      +ux_temp*uz_temp*cx[pos]*cz[pos]
										      +uy_temp*uz_temp*cy[pos]*cz[pos])))
		+kconst*(wxx[pos]*gradx_temp*gradx_temp+wyy[pos]*grady_temp*grady_temp
		        +wzz[pos]*gradz_temp*gradz_temp+wxy[pos]*gradx_temp*grady_temp
		        +wyz[pos]*grady_temp*gradz_temp+wzx[pos]*gradx_temp*gradz_temp);

   		double phasepos=weights[pos]*(chemical_pot+phase_press*(cx[pos]*ux_temp+cy[pos]*uy_temp+cz[pos]*uz_temp)
						+1.5*phase_press*((cx[pos]*cx[pos]-1.0/3.0)*ux_temp*ux_temp
                                        +(cy[pos]*cy[pos]-1.0/3.0)*uy_temp*uy_temp
                                        +(cz[pos]*cz[pos]-1.0/3.0)*uz_temp*uz_temp
                                        +2.0*(ux_temp*uy_temp*cx[pos]*cy[pos]
                                             +ux_temp*uz_temp*cx[pos]*cz[pos]
                                             +uy_temp*uz_temp*cy[pos]*cz[pos])));


        int neg=unknown[compliment[iCount]];
        double eqneg=weights[neg]*(pressure_bulk+rho_press*(cx[neg]*ux_temp+cy[neg]*uy_temp+cz[neg]*uz_temp)
						+1.5*rho_press*((cx[neg]*cx[neg]-1.0/3.0)*ux_temp*ux_temp
                                         +(cy[neg]*cy[neg]-1.0/3.0)*uy_temp*uy_temp
                                         +(cz[neg]*cz[neg]-1.0/3.0)*uz_temp*uz_temp
										 +2.0*(ux_temp*uy_temp*cx[neg]*cy[neg]
										      +ux_temp*uz_temp*cx[neg]*cz[neg]
										      +uy_temp*uz_temp*cy[neg]*cz[neg])))
		+kconst*(wxx[neg]*gradx_temp*gradx_temp+wyy[neg]*grady_temp*grady_temp
		        +wzz[neg]*gradz_temp*gradz_temp+wxy[neg]*gradx_temp*grady_temp
		        +wyz[neg]*grady_temp*gradz_temp+wzx[neg]*gradx_temp*gradz_temp);

        double phaseneg=weights[neg]*(chemical_pot+phase_press*(cx[neg]*ux_temp+cy[neg]*uy_temp+cz[neg]*uz_temp)
						+1.5*phase_press*((cx[neg]*cx[neg]-1.0/3.0)*ux_temp*ux_temp
                                        +(cy[neg]*cy[neg]-1.0/3.0)*uy_temp*uy_temp
                                        +(cz[neg]*cz[neg]-1.0/3.0)*uz_temp*uz_temp
                                        +2.0*(ux_temp*uy_temp*cx[neg]*cy[neg]
                                             +ux_temp*uz_temp*cx[neg]*cz[neg]
                                             +uy_temp*uz_temp*cy[neg]*cz[neg])));

        //Perform BB for the Non-equilibrium parts
        lattice->f[counter*NPOP+iCount]=lattice->f[counter*NPOP+compliment[iCount]]+(eqpos-eqneg);

        //Do we need to perform the same for the phase population?
        lattice->g[counter*NPOP+iCount]=lattice->g[counter*NPOP+compliment[iCount]]+(phasepos-phaseneg);

    }

    //Calculate the new densities and velocities
    double phase_temp=0.0, dense_temp=0.0;

    ux_temp=0.0;
    uy_temp=0.0;
    uz_temp=0.0;

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

    double diff[3];

    diff[0]=(lattice->rho[counter]*lattice->ux[counter]-dense_temp*ux_temp)/double(diagonal.size());
    diff[1]=(lattice->rho[counter]*lattice->uy[counter]-dense_temp*uy_temp)/double(diagonal.size());
    diff[2]=(lattice->rho[counter]*lattice->uz[counter]-dense_temp*uz_temp)/double(diagonal.size());

    for(int iPop=0;iPop<diagonal.size();iPop++)
    {
       lattice->f[counter*NPOP+diagonal[iPop]] +=
               (!normx)*cx[diagonal[iPop]]*diff[0]
               +(!normy)*cy[diagonal[iPop]]*diff[1]
               +(!normz)*cz[diagonal[iPop]]*diff[2];
    }


	//double phase_square=phase_press*phase_press;
	//double pressure_bulk=dense_press/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_press*laplace_temp;
	//double chemical_pot=gammaconst*(aconst*(-phase_press+phase_press*phase_press*phase_press)-kconst*laplace_temp);

    double sum=0.0;
	double sum_phase=0.0;

	double feq[NPOP], geq[NPOP];
	for (int k=1; k<NPOP; k++)
	{
			feq[k]=weights[k]*(pressure_bulk+rho_press*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)
							+1.5*rho_press*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+(cz[k]*cz[k]-1.0/3.0)*uz_temp*uz_temp
												 +2.0*(ux_temp*uy_temp*cx[k]*cy[k]+ux_temp*uz_temp*cx[k]*cz[k]+uy_temp*uz_temp*cy[k]*cz[k])))
			+kconst*(wxx[k]*gradx_temp*gradx_temp+wyy[k]*grady_temp*grady_temp+wzz[k]*gradz_temp*gradz_temp
						 +wxy[k]*gradx_temp*grady_temp+wyz[k]*grady_temp*gradz_temp+wzx[k]*gradx_temp*gradz_temp);
			geq[k]=weights[k]*(chemical_pot+phase_press*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)
							+1.5*phase_press*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+(cz[k]*cz[k]-1.0/3.0)*uz_temp*uz_temp
												 +2.0*(ux_temp*uy_temp*cx[k]*cy[k]+ux_temp*uz_temp*cx[k]*cz[k]+uy_temp*uz_temp*cy[k]*cz[k])));

			sum+=feq[k];
			sum_phase+=geq[k];
	}

	feq[0]=rho_press-sum;
	geq[0]=phase_press-sum_phase;

	//Omega calculation
	double tau_temp=tau_gas+(phase_temp+1.0)/2.0*(tau_liq-tau_gas);
	double omega_temp=1.0/tau_temp;

	//Collision procedure
	for (int k=0; k<NPOP; k++)
	{
		//cout<<"Increment="<<-omega*(lattice->f[counter*NPOP+k]-feq[k])<<"\n";
		lattice->f[counter*NPOP+k]+=-omega_temp*(lattice->f[counter*NPOP+k]-feq[k]);
		lattice->g[counter*NPOP+k]+=-omega_temp*(lattice->g[counter*NPOP+k]-geq[k]);
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



void DynamicsZouHePressure::init(int iX,int iY, int iZ)
{
    //cout<<"We are here \n";

    Lattice * lattice=this->solver->getLattice();

    int NX=this->solver->geom->getNX();
    int NY=this->solver->geom->getNY();

    int counter=iZ*NX*NY+iY*NX+iX;

   	int zbegin=lattice->getZbegin();
	int zend=lattice->getZend();

	//Calculation of gradients and so on
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

  	double sum=0.0;
	double sum_phase=0.0;
	double phase_square=phase_press*phase_press;
	double pressure_bulk=rho_press/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_press*laplace_temp;
	double chemical_pot=gammaconst*(aconst*(-phase_press+phase_press*phase_press*phase_press)-kconst*laplace_temp);

 	double feq[NPOP], geq[NPOP];
	for (int k=1; k<NPOP; k++)
	{
			feq[k]=weights[k]*(pressure_bulk+rho_press*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)
							+1.5*rho_press*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+(cz[k]*cz[k]-1.0/3.0)*uz_temp*uz_temp
												 +2.0*(ux_temp*uy_temp*cx[k]*cy[k]+ux_temp*uz_temp*cx[k]*cz[k]+uy_temp*uz_temp*cy[k]*cz[k])))
			+kconst*(wxx[k]*gradx_temp*gradx_temp+wyy[k]*grady_temp*grady_temp+wzz[k]*gradz_temp*gradz_temp
						 +wxy[k]*gradx_temp*grady_temp+wyz[k]*grady_temp*gradz_temp+wzx[k]*gradx_temp*gradz_temp);
			geq[k]=weights[k]*(chemical_pot+phase_press*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)
							+1.5*phase_press*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+(cz[k]*cz[k]-1.0/3.0)*uz_temp*uz_temp
												 +2.0*(ux_temp*uy_temp*cx[k]*cy[k]+ux_temp*uz_temp*cx[k]*cz[k]+uy_temp*uz_temp*cy[k]*cz[k])));

			sum+=feq[k];
			sum_phase+=geq[k];
	}

	feq[0]=rho_press-sum;
	geq[0]=phase_press-sum_phase;

    for(int iPop=0; iPop<NPOP;iPop++)
    {
        lattice->f[counter*NPOP+iPop]=feq[iPop];
        lattice->g[counter*NPOP+iPop]=geq[iPop];
    }

}

