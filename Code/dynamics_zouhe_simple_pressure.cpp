#include "dynamics_zouhe_simple_pressure.h"
#include "solver.h"

DynamicsZouHeSimplePressure::DynamicsZouHeSimplePressure(Solver * _solver,ParamsList _params_list):
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
            if ( ! ((cx[iCoor]==normx) && (cy[iCoor]==normy) && (cz[iCoor]==normz)) )
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

void DynamicsZouHeSimplePressure::updateFields(int iX,int iY,int iZ)
{
  	int NX=this->solver->geom->getNX();
	int NY=this->solver->geom->getNY();
    Lattice * lattice=this->solver->getLattice();


    int counter=iZ*NX*NY+iY*NX+iX;
    lattice->rho[counter]=rho_press;
    lattice->phase[counter]=phase_press;

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

    double norma=sqrt(normx*normx+normy*normy+normz*normz);

    //Commented here
    //double velocity_perp=((-2.0*denseminus-denseneutral+rho_press)/rho_press)*
    //                        (double(normx)+double(normy)+double(normz))/norma;

    double velocity_perp=((-2.0*denseminus-denseneutral+rho_press)/rho_press);

    double velocity_tang=0.0;

    lattice->ux[counter]=velocity_perp*double(normx)/norma;
    lattice->uy[counter]=velocity_perp*double(normy)/norma;
    lattice->uz[counter]=velocity_perp*double(normz)/norma;

}

void DynamicsZouHeSimplePressure::collide_stream(int iX,int iY,int iZ)
{
    int NX=this->solver->geom->getNX();
    int NY=this->solver->geom->getNY();
    Lattice * lattice=this->solver->getLattice();
   	int zbegin=lattice->getZbegin();
	int zend=lattice->getZend();
    int counter=iZ*NX*NY+iY*NX+iX;


    //Testing scheme
    double pop_test[19],pop_after[19];
    for(int iPop=0;iPop<19;iPop++)
        pop_test[iPop]=pop_after[iPop]=lattice->f[counter*NPOP+iPop];



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

    double norma=sqrt(normx*normx+normy*normy+normz*normz);

    //Commented here
    //double velocity_perp=((-2.0*denseminus-denseneutral+rho_press)/rho_press)*
    //                        (double(normx)+double(normy)+double(normz))/norma;
    double velocity_perp=((-2.0*denseminus-denseneutral+rho_press)/rho_press);
    if (iZ==0)
    {
        //cout<<"velocity_perp="<<velocity_perp<<"\n";
    }

    double velocity_tang=0.0;

    //if (normz==1)
    //    cout<<"We see negative as well\n";
    //Transform velocities to real velocities
    //cout<<"norma="<<norma<<"\n";
    double ux_temp=velocity_perp*normx/norma;
    double uy_temp=velocity_perp*normy/norma;
    double uz_temp=velocity_perp*normz/norma;


    double ux_test=ux_temp;
    double uy_test=uy_temp;
    double uz_test=uz_temp;
    //cout<<"uz_temp"<<uz_temp<<"\n";

    for(int iCount=0;iCount<unknown.size();iCount++)
    {
        int pos=unknown[iCount];

        double eqpos=weights[pos]*(rho_press/3.0+rho_press*(cx[pos]*ux_temp+cy[pos]*uy_temp+cz[pos]*uz_temp)
						+1.5*rho_press*((cx[pos]*cx[pos]-1.0/3.0)*ux_temp*ux_temp
                                         +(cy[pos]*cy[pos]-1.0/3.0)*uy_temp*uy_temp
                                         +(cz[pos]*cz[pos]-1.0/3.0)*uz_temp*uz_temp
										 +2.0*(ux_temp*uy_temp*cx[pos]*cy[pos]
										      +ux_temp*uz_temp*cx[pos]*cz[pos]
										      +uy_temp*uz_temp*cy[pos]*cz[pos])));

   		double phasepos=weights[pos]*(phase_press/3.0+phase_press*(cx[pos]*ux_temp+cy[pos]*uy_temp+cz[pos]*uz_temp)
						+1.5*phase_press*((cx[pos]*cx[pos]-1.0/3.0)*ux_temp*ux_temp
                                        +(cy[pos]*cy[pos]-1.0/3.0)*uy_temp*uy_temp
                                        +(cz[pos]*cz[pos]-1.0/3.0)*uz_temp*uz_temp
                                        +2.0*(ux_temp*uy_temp*cx[pos]*cy[pos]
                                             +ux_temp*uz_temp*cx[pos]*cz[pos]
                                             +uy_temp*uz_temp*cy[pos]*cz[pos])));

        int neg=compliment[unknown[iCount]];
        double eqneg=weights[neg]*(rho_press/3.0+rho_press*(cx[neg]*ux_temp+cy[neg]*uy_temp+cz[neg]*uz_temp)
						+1.5*rho_press*((cx[neg]*cx[neg]-1.0/3.0)*ux_temp*ux_temp
                                         +(cy[neg]*cy[neg]-1.0/3.0)*uy_temp*uy_temp
                                         +(cz[neg]*cz[neg]-1.0/3.0)*uz_temp*uz_temp
										 +2.0*(ux_temp*uy_temp*cx[neg]*cy[neg]
										      +ux_temp*uz_temp*cx[neg]*cz[neg]
										      +uy_temp*uz_temp*cy[neg]*cz[neg])));

        double phaseneg=weights[neg]*(phase_press/3.0+phase_press*(cx[neg]*ux_temp+cy[neg]*uy_temp+cz[neg]*uz_temp)
						+1.5*phase_press*((cx[neg]*cx[neg]-1.0/3.0)*ux_temp*ux_temp
                                        +(cy[neg]*cy[neg]-1.0/3.0)*uy_temp*uy_temp
                                        +(cz[neg]*cz[neg]-1.0/3.0)*uz_temp*uz_temp
                                        +2.0*(ux_temp*uy_temp*cx[neg]*cy[neg]
                                             +ux_temp*uz_temp*cx[neg]*cz[neg]
                                             +uy_temp*uz_temp*cy[neg]*cz[neg])));

        //Perform BB for the Non-equilibrium parts
        lattice->f[counter*NPOP+pos]=lattice->f[counter*NPOP+neg]+(eqpos-eqneg);

        //cout<<pos<<"\n";
        if (pos==16)
        {
            //cout<<"ux="<<ux_temp<<" uy="<<uy_temp<<" uz="<<uz_temp<<"\n";
            //cout<<"Diff="<<eqpos-eqneg<<"should be "<<2.0/12.0*rho_press*uz_temp<<"\n";
        }

        //Do we need to perform the same for the phase population?
        lattice->g[counter*NPOP+pos]=lattice->g[counter*NPOP+neg]+(phasepos-phaseneg);

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

    diff[0]=-dense_temp*ux_temp/double(diagonal.size());
    diff[1]=-dense_temp*uy_temp/double(diagonal.size());
    diff[2]=-dense_temp*uz_temp/double(diagonal.size());

    if (iZ==0)
    {
//        cout<<"Compare this stuff with another\n";
//        cout<<"X-component="<<diff[0]<<"compare with calculated "<<
//              0.5*(lattice->f[counter*NPOP+1]+lattice->f[counter*NPOP+7]+lattice->f[counter*NPOP+8]
//              -(lattice->f[counter*NPOP+2]+lattice->f[counter*NPOP+11]+lattice->f[counter*NPOP+12]))<<"\n";
//        cout<<"Y-component="<<diff[2]<<"compare with calculated "<<
//              0.5*(lattice->f[counter*NPOP+3]+lattice->f[counter*NPOP+7]+lattice->f[counter*NPOP+11]
//              -(lattice->f[counter*NPOP+4]+lattice->f[counter*NPOP+8]+lattice->f[counter*NPOP+12]))<<"\n";

    }

    //cout<<"!normx "<<!normx<<" !normy "<<!normy<<" !normz "<<!normz<<"\n";

    for(int iPop=0;iPop<diagonal.size();iPop++)
    {
       lattice->f[counter*NPOP+diagonal[iPop]] +=
               (!normx)*cx[diagonal[iPop]]*diff[0]
               +(!normy)*cy[diagonal[iPop]]*diff[1]
               +(!normz)*cz[diagonal[iPop]]*diff[2];
    }


    for(int counter=0;counter<unknown.size();counter++)
    {
        int iPop=unknown[counter];
        int tx=cx[iPop]-(cx[iPop]*normx+cy[iPop]*normy+cz[iPop]*normz)*normx;
        int ty=cy[iPop]-(cx[iPop]*normx+cy[iPop]*normy+cz[iPop]*normz)*normy;
        int tz=cz[iPop]-(cx[iPop]*normx+cy[iPop]*normy+cz[iPop]*normz)*normz;

        double newpop=pop_test[iPop];
        for (int jPop=1;jPop<19;jPop++)
        {
            newpop+=0.5*pop_test[jPop]*(tx*cx[jPop]+ty*cy[jPop]+tz*cz[jPop])*
                    (1.0-sqrt(cx[jPop]*cx[jPop]*normx*normx+cy[jPop]*cy[jPop]*normy*normy+cz[jPop]*cz[jPop]*normz*normz));
        }
        newpop-=rho_press/6.0*(cx[iPop]*ux_test+cy[iPop]*uy_test+cz[iPop]*uz_test)
                +rho_press/3.0*(tx*ux_test+ty*uy_test+tz*uz_test);
        //cout<<newpop<<"compare this stuff with "<<lattice->f[counter*NPOP+iPop]<<"\n";
        pop_after[iPop]=newpop;
    }

//    //change to newpop
//    for(int iPop=0;iPop<19;iPop++)
//    {
//        lattice->f[counter*NPOP+iPop]=pop_after[iPop];
//    }

    //To calculate the previous stuff
    ux_temp=velocity_perp*normx/norma;
    uy_temp=velocity_perp*normy/norma;
    uz_temp=velocity_perp*normz/norma;


    //Equilibriun function construction
    double sum=0.0;
	double sum_phase=0.0;

	double feq[NPOP], geq[NPOP];
	for (int k=1; k<NPOP; k++)
	{
			feq[k]=weights[k]*(rho_press/3.0+rho_press*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)
							+1.5*rho_press*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+(cz[k]*cz[k]-1.0/3.0)*uz_temp*uz_temp
												 +2.0*(ux_temp*uy_temp*cx[k]*cy[k]+ux_temp*uz_temp*cx[k]*cz[k]+uy_temp*uz_temp*cy[k]*cz[k])));
			geq[k]=weights[k]*(phase_press/3.0+phase_press*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)
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



void DynamicsZouHeSimplePressure::init(int iX,int iY, int iZ)
{
    //cout<<"We are here \n";

    Lattice * lattice=this->solver->getLattice();

    int NX=this->solver->geom->getNX();
    int NY=this->solver->geom->getNY();

    int counter=iZ*NX*NY+iY*NX+iX;

   	int zbegin=lattice->getZbegin();
	int zend=lattice->getZend();

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
			feq[k]=weights[k]*(rho_press/3.0+rho_press*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)
							+1.5*rho_press*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+(cz[k]*cz[k]-1.0/3.0)*uz_temp*uz_temp
												 +2.0*(ux_temp*uy_temp*cx[k]*cy[k]+ux_temp*uz_temp*cx[k]*cz[k]+uy_temp*uz_temp*cy[k]*cz[k])));
			geq[k]=weights[k]*(phase_press/3.0+phase_press*(cx[k]*ux_temp+cy[k]*uy_temp+cz[k]*uz_temp)
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

