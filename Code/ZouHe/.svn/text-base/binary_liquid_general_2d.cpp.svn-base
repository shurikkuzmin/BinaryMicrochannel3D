#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>

//Domain size
const int NY=9;
const int NX=101;

//Time steps
const int N=1000;
const int NOUTPUT=100;

//Fields and populations
double f[NX][NY][9], f2[NX][NY][9], g[NX][NY][9], g2[NX][NY][9];
double rho[NX][NY],ux[NX][NY],uy[NX][NY],phase[NX][NY];

//Pressure boundary conditions
double rho_inlet=1.003;
double rho_outlet=1.0;
double phase_inlet=1.0;
double phase_outlet=1.0;

//Binary-liquid parameters
double aconst=0.04;
double kconst=0.04;
double gammaconst=1.0;

//BGK relaxation parameter
double omega=1.0;

//Magic Irina's parameters
double omegaginzburg=8.0*(2.0-omega)/(8.0-omega);
double omegamat[]={1.0,1.0,1.0,omega,omega,omega,1.0,omegaginzburg,omegaginzburg};

//Underlying lattice parameters
double weights[]={4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};
int cx[]={0,1,0,-1,0,1,-1,-1,1};
int cy[]={0,0,1,0,-1,1,1,-1,-1};
double gmat[]={1.0,-2.0,-2.0,-2.0,-2.0,4.0,4.0,4.0,4.0};
int compliment[]={0,3,4,1,2,7,8,5,6};
float wxx[] = {0.0, 1.0/3.0, -1.0/6.0, 1.0/3.0, -1.0/6.0, -1.0/24.0, -1.0/24.0, -1.0/24.0, -1.0/24.0};
float wyy[] = {0.0, -1.0/6.0, 1.0/3.0, -1.0/6.0, 1.0/3.0, -1.0/24.0, -1.0/24.0, -1.0/24.0, -1.0/24.0};
float wxy[] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0/4.0, -1.0/4.0, 1.0/4.0, -1.0/4.0};

float gradstencilx[9]={0.0,4.0/12.0,0.0,-4.0/12.0,0.0,
			          1.0/12.0,-1.0/12.0,-1.0/12.0,1.0/12.0};

float gradstencily[9]={0.0,0.0,4.0/12.0,0.0,-4.0/12.0,
			          1.0/12.0,1.0/12.0,-1.0/12.0,-1.0/12.0};

float laplacestencil[9]={-20.0/6.0,4.0/6.0,4.0/6.0,4.0/6.0,4.0/6.0,
					   1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0};

//Matrix M
double M[9][9];


void matrix_init()
{
	for (int iCoor=0;iCoor<9;iCoor++)
	{
		M[0][iCoor]=1.0;
		M[1][iCoor]=cx[iCoor]*sqrt(3.0);
		M[2][iCoor]=cy[iCoor]*sqrt(3.0);
		M[3][iCoor]=(cx[iCoor]*cx[iCoor]-1.0/3.0)*3.0/sqrt(2.0);
		M[4][iCoor]=cx[iCoor]*cy[iCoor]*3.0;
		M[5][iCoor]=(cy[iCoor]*cy[iCoor]-1.0/3.0)*3.0/sqrt(2.0);
		M[6][iCoor]=gmat[iCoor]/2.0;
		M[7][iCoor]=gmat[iCoor]*cx[iCoor]*sqrt(1.5)/2.0;
		M[8][iCoor]=gmat[iCoor]*cy[iCoor]*sqrt(1.5)/2.0;
	}
}

void writedensity(std::string const & fName)
{
	std::string fullName = "./tmp/" + fName+ ".dat";
	std::ofstream fout(fullName.c_str());
	fout.precision(10);

	for (int iY=NY-1; iY>=0; --iY)
	{
		for (int iX=0; iX<NX; ++iX)
			fout<<rho[iX][iY]<<" ";
		fout<<"\n";
	}

}

void writephase(std::string const & fName)
{
	std::string fullName = "./tmp/" + fName+ ".dat";
	std::ofstream fout(fullName.c_str());
	fout.precision(10);

	for (int iY=NY-1; iY>=0; --iY)
	{
		for (int iX=0; iX<NX; ++iX)
			fout<<phase[iX][iY]<<" ";
		fout<<"\n";
	}

}


void writevelocity(std::string const & fName)
{
	std::string fullName = "./tmp/" + fName+ ".dat";
	std::ofstream fout(fullName.c_str());
	fout.precision(10);

	for (int iY=NY-1; iY>=0; --iY)
	{
		for (int iX=0; iX<NX; ++iX)
			fout<<ux[iX][iY]<<" ";
		fout<<"\n";
	}

}

void init()
{

	//Phase initialization prior any equilibrium functions calculations
    for(int iX=0;iX<NX;iX++)
		for(int iY=0; iY<NY; iY++)
		{
			if ( (iX>=(NX-1)/4) && (iX<=3*(NX-1)/4) && (iY>=2) && (iY<=NY-3) )
            {
                phase[iX][iY]=-1.0;
            }
            else
            {
                phase[iX][iY]=1.0;
            }

			if (iX==0)
			{
				phase[iX][iY]=phase_inlet;
			}

			if (iX==NX-1)
			{
				phase[iX][iY]=phase_outlet;
			}

			if ((iY==0)||(iY==NY-1))
				phase[iX][iY]=0.5;

		}

	//Bulk nodes initialization
	for(int iX=1;iX<NX-1;iX++)
		for(int iY=1;iY<NY-1;iY++)
		{


			double laplace_temp=0.0;
			double gradx_temp=0.0;
			double grady_temp=0.0;
			for(int k=0;k<9;k++)
			{
				int iX2=(iX+cx[k]+NX) % NX;
				int iY2=(iY+cy[k]+NY) % NY;
				laplace_temp+=laplacestencil[k]*phase[iX2][iY2];
				gradx_temp+=gradstencilx[k]*phase[iX2][iY2];
				grady_temp+=gradstencily[k]*phase[iX2][iY2];
			}


			//Initialization of the macroscopic fields
			rho[iX][iY]=1.0;
			ux[iX][iY]=0.0;
			uy[iX][iY]=0.0;


			double phase_temp=phase[iX][iY];
			double dense_temp=rho[iX][iY];
			double ux_temp=ux[iX][iY];
			double uy_temp=uy[iX][iY];

			double feq;
			double geq;
			double sum=0.0;
			double sum_phase=0.0;
			double phase_square=phase_temp*phase_temp;
			double pressure_bulk=dense_temp/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_temp*laplace_temp;
			double chemical_pot=gammaconst*(aconst*(-phase_temp+phase_temp*phase_temp*phase_temp)-kconst*laplace_temp);

			for (int k=1; k<9; k++)
			{
				feq=weights[k]*(3.0*pressure_bulk+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+2.0*ux_temp*uy_temp*cx[k]*cy[k]))
					+kconst*(wxx[k]*gradx_temp*gradx_temp+wyy[k]*grady_temp*grady_temp+wxy[k]*gradx_temp*grady_temp);
				geq=weights[k]*(3.0*chemical_pot+3.0*phase_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*phase_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp
												 +2.0*ux_temp*uy_temp*cx[k]*cy[k]));
				sum+=feq;
				sum_phase+=geq;

				f[iX][iY][k]=feq;
				g[iX][iY][k]=geq;
			}
			f[iX][iY][0]=dense_temp-sum;
			g[iX][iY][0]=phase_temp-sum_phase;

		}

	//Pressure boundary nodes initialization
	for(int iY=1;iY<NY-1;iY++)
	{
		rho[0][iY]=rho_inlet;
		ux[0][iY]=0.0;
		uy[0][iY]=0.0;

		double fluxx=3.0*ux[0][iY];
		double fluxy=3.0*uy[0][iY];

		double qxx=4.5*ux[0][iY]*ux[0][iY];
		double qxy=9.0*ux[0][iY]*uy[0][iY];
		double qyy=4.5*uy[0][iY]*uy[0][iY];

		for(int iPop=0;iPop<9;iPop++)
		{
			f[0][iY][iPop]=weights[iPop]*rho[0][iY]*(1.0+fluxx*cx[iPop]+fluxy*cy[iPop]+
							qxx*(cx[iPop]*cx[iPop]-1.0/3.0)+qxy*cx[iPop]*cy[iPop]+qyy*(cy[iPop]*cy[iPop]-1.0/3.0));
			g[0][iY][iPop]=weights[iPop]*phase[0][iY]*(1.0+fluxx*cx[iPop]+fluxy*cy[iPop]+
										   qxx*(cx[iPop]*cx[iPop]-1.0/3.0)+qxy*cx[iPop]*cy[iPop]+qyy*(cy[iPop]*cy[iPop]-1.0/3.0));
		}


   		rho[NX-1][iY]=rho_outlet;
		ux[NX-1][iY]=0.0;
		uy[NX-1][iY]=0.0;

		fluxx=3.0*ux[NX-1][iY];
		fluxy=3.0*uy[NX-1][iY];

		qxx=4.5*ux[NX-1][iY]*ux[NX-1][iY];
		qxy=9.0*ux[NX-1][iY]*uy[NX-1][iY];
		qyy=4.5*uy[NX-1][iY]*uy[NX-1][iY];

		for(int iPop=0;iPop<9;iPop++)
		{
			f[NX-1][iY][iPop]=weights[iPop]*rho[NX-1][iY]*(1.0+fluxx*cx[iPop]+fluxy*cy[iPop]+
							qxx*(cx[iPop]*cx[iPop]-1.0/3.0)+qxy*cx[iPop]*cy[iPop]+qyy*(cy[iPop]*cy[iPop]-1.0/3.0));
			g[NX-1][iY][iPop]=weights[iPop]*phase[NX-1][iY]*(1.0+fluxx*cx[iPop]+fluxy*cy[iPop]+
														   qxx*(cx[iPop]*cx[iPop]-1.0/3.0)+qxy*cx[iPop]*cy[iPop]+qyy*(cy[iPop]*cy[iPop]-1.0/3.0));

		}
	}

	//BB nodes initialization
    for(int iX=0;iX<NX;iX++)
	{
		rho[iX][0]=1.0;
		phase[iX][0]=0.5;
		ux[iX][0]=0.0;
		uy[iX][0]=0.0;

		double fluxx=3.0*ux[iX][0];
		double fluxy=3.0*uy[iX][0];

		double qxx=4.5*ux[iX][0]*ux[iX][0];
		double qxy=9.0*ux[iX][0]*uy[iX][0];
		double qyy=4.5*uy[iX][0]*uy[iX][0];


		for(int iPop=0;iPop<9;iPop++)
		{
			f[iX][0][iPop]=weights[iPop]*rho[iX][0]*(1.0+fluxx*cx[iPop]+fluxy*cy[iPop]+
							qxx*(cx[iPop]*cx[iPop]-1.0/3.0)+qxy*cx[iPop]*cy[iPop]+qyy*(cy[iPop]*cy[iPop]-1.0/3.0));
			g[iX][0][iPop]=weights[iPop]*phase[iX][0]*(1.0+fluxx*cx[iPop]+fluxy*cy[iPop]+
														   qxx*(cx[iPop]*cx[iPop]-1.0/3.0)+qxy*cx[iPop]*cy[iPop]+qyy*(cy[iPop]*cy[iPop]-1.0/3.0));

		}

		rho[iX][NY-1]=1.0;
		phase[iX][NY-1]=0.5;
        ux[iX][NY-1]=0.0;
        uy[iX][NY-1]=0.0;

		fluxx=3.0*ux[iX][NY-1];
		fluxy=3.0*uy[iX][NY-1];

		qxx=4.5*ux[iX][NY-1]*ux[iX][NY-1];
		qxy=9.0*ux[iX][NY-1]*uy[iX][NY-1];
		qyy=4.5*uy[iX][NY-1]*uy[iX][NY-1];


		for(int iPop=0;iPop<9;iPop++)
		{
			f[iX][NY-1][iPop]=weights[iPop]*rho[iX][NY-1]*(1.0+fluxx*cx[iPop]+fluxy*cy[iPop]+
							qxx*(cx[iPop]*cx[iPop]-1.0/3.0)+qxy*cx[iPop]*cy[iPop]+qyy*(cy[iPop]*cy[iPop]-1.0/3.0));
			g[iX][NY-1][iPop]=weights[iPop]*phase[iX][NY-1]*(1.0+fluxx*cx[iPop]+fluxy*cy[iPop]+
														   qxx*(cx[iPop]*cx[iPop]-1.0/3.0)+qxy*cx[iPop]*cy[iPop]+qyy*(cy[iPop]*cy[iPop]-1.0/3.0));

		}


	}

}

void general_pressure(int xCoor, int normx, int normy,double rho_pressure)
{
    std::vector<int> unknown;
    std::vector<int> known;
    std::vector<int> positive;
    std::vector<int> negative;
    std::vector<int> neutral;
    std::vector<int> diagonal;

    for(int iCoor=0;iCoor<9;iCoor++)
    {
        //Determination of the projection
        int dot_product=cx[iCoor]*normx+cy[iCoor]*normy;

        if (dot_product>0)
        {
            unknown.push_back(iCoor);
            positive.push_back(iCoor);
            if (! ((cx[iCoor]==normx)&&(cy[iCoor]==normy)) )
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

//    for(int iCoor=0;iCoor<neutral.size();iCoor++)
//    {
//        std::cout<<"Neutral populations="<<neutral[iCoor]<<"\n";
//    }
//    std::cout<<"\n";

    //Update macroscopic fields
    for (int iY=1;iY<NY-1;iY++)
    {
        rho[xCoor][iY]=rho_pressure;

        double denseminus=0.0;
        double denseneutral=0.0;

        for(int iCount=0;iCount<negative.size();iCount++)
            denseminus+=f[xCoor][iY][negative[iCount]];

        for(int iCount=0;iCount<neutral.size();iCount++)
            denseneutral+=f[xCoor][iY][neutral[iCount]];

        double norma=sqrt(normx*normx+normy*normy);

        double velocity_perp=(-2.0*denseminus-denseneutral+rho_pressure)/rho_pressure;

        ux[xCoor][iY]=velocity_perp*double(normx)/norma;
        uy[xCoor][iY]=velocity_perp*double(normy)/norma;

        double rho_temp=rho[xCoor][iY];
        double ux_temp=ux[xCoor][iY];
        double uy_temp=uy[xCoor][iY];

        //Perform non-equilibrium BB for the unknown populations
        for(int iCount=0;iCount<unknown.size();iCount++)
        {
            int pos=unknown[iCount];


            double eqpos=weights[pos]*(rho_temp + 3.0*rho_temp*(cx[pos]*ux_temp+cy[pos]*uy_temp)
                            +4.5*rho_temp*((cx[pos]*cx[pos]-1.0/3.0)*ux_temp*ux_temp
                                             +(cy[pos]*cy[pos]-1.0/3.0)*uy_temp*uy_temp
                                             +2.0*ux_temp*uy_temp*cx[pos]*cy[pos]));

            int neg=compliment[unknown[iCount]];
            double eqneg=weights[neg]*(rho_temp + 3.0*rho_temp*(cx[neg]*ux_temp+cy[neg]*uy_temp)
                            +4.5*rho_temp*((cx[neg]*cx[neg]-1.0/3.0)*ux_temp*ux_temp
                                             +(cy[neg]*cy[neg]-1.0/3.0)*uy_temp*uy_temp
                                             +2.0*ux_temp*uy_temp*cx[neg]*cy[neg]));
            //Perform BB for the Non-equilibrium parts
            f[xCoor][iY][pos]=f[xCoor][iY][neg]+(eqpos-eqneg);
        }

        //Calculate the excess of momenta
        double momx=0.0;
        double momy=0.0;

        for(int iPop=0;iPop<9;iPop++)
        {
			momx+=f[xCoor][iY][iPop]*cx[iPop];
			momy+=f[xCoor][iY][iPop]*cy[iPop];
		}

        momx=(rho_temp*ux_temp-momx)/diagonal.size();
        momy=(rho_temp*uy_temp-momy)/diagonal.size();
        //std::cout<<f[xCoor][iY][2]-f[xCoor][iY][4]<<" compare with "<<2*momy<<"\n";


        //Add the correction to the diagonal populations
        for(int iPop=0;iPop<diagonal.size();iPop++)
        {
            f[xCoor][iY][diagonal[iPop]] +=
               (!normx)*cx[diagonal[iPop]]*momx
               +(!normy)*cy[diagonal[iPop]]*momy;
        }

		double feqeq[9];
		for(int iPop=0;iPop<9; iPop++)
		{
// 					feqeq[iPop]=weights[iPop]*dense*(1.0+3.0*cx[iPop]*uxeq+3.0*cy[iPop]*uyeq+
// 					4.5*(cx[iPop]*cx[iPop]-1.0/3.0)*uxeq*uxeq+9.0*cx[iPop]*cy[iPop]*uxeq*uyeq+4.5*(cy[iPop]*cy[iPop]-1.0/3.0)*uyeq*uyeq+
// 					81.0/4.0*(cx[iPop]*cx[iPop]*cy[iPop]*cy[iPop]-1.0/3.0*cx[iPop]*cx[iPop]-1.0/3.0*cy[iPop]*cy[iPop]+1.0/9.0)*uxeq*uxeq*uyeq*uyeq+
// 					27.0/2.0*(cx[iPop]*cx[iPop]*cy[iPop]-1.0/3.0*cy[iPop])*uxeq*uxeq*uyeq+
// 					27.0/2.0*(cy[iPop]*cy[iPop]*cx[iPop]-1.0/3.0*cx[iPop])*uyeq*uyeq*uxeq);

			feqeq[iPop]=weights[iPop]*rho_temp*(1.0+3.0*cx[iPop]*ux_temp+3.0*cy[iPop]*uy_temp+
                        4.5*(cx[iPop]*cx[iPop]-1.0/3.0)*ux_temp*ux_temp
                        +9.0*cx[iPop]*cy[iPop]*ux_temp*uy_temp
                        +4.5*(cy[iPop]*cy[iPop]-1.0/3.0)*uy_temp*uy_temp);

			f2[xCoor][iY][iPop]=f[xCoor][iY][iPop]*(1.0-omega)+omega*feqeq[iPop];
        }

		//For the time being we initialize the phase populations through the equilibrium function

		double phase_temp;
		int xCoor2;
		if(xCoor==0)
		{
			phase[xCoor][iY]=phase_inlet;
			phase_temp=phase_inlet;
            xCoor2=1;

            //Outflow BC
            //phase[xCoor][iY]=phase[xCoor2][iY];
			//phase_temp=phase[xCoor][iY];

		}
		else
		{
			phase[xCoor][iY]=phase_outlet;
			phase_temp=phase_outlet;
            xCoor2=NX-2;

            //Outflow BC
            //phase[xCoor][iY]=phase[xCoor2][iY];
			//phase_temp=phase[xCoor][iY];

		}
		double geqeq[9];

		for(int iPop=0;iPop<9; iPop++)
		{

			//Full Equilibrium
//			geqeq[iPop]=weights[iPop]*phase_temp*(1.0+3.0*cx[iPop]*ux_temp+3.0*cy[iPop]*uy_temp
//                                       +4.5*(cx[iPop]*cx[iPop]-1.0/3.0)*ux_temp*ux_temp
//									   +9.0*cx[iPop]*cy[iPop]*ux_temp*uy_temp
//									   +4.5*(cy[iPop]*cy[iPop]-1.0/3.0)*uy_temp*uy_temp);

            //No velocity equilibrium
            geqeq[iPop]=weights[iPop]*phase_temp;

            //Equilibrium populations initialization - tends to be unstable
			g2[xCoor][iY][iPop]=geqeq[iPop];

			//OutFlow conditions
			//g[xCoor][iY][iPop]=g[xCoor2][iY][iPop];
			//g2[xCoor][iY][iPop]=g[xCoor][iY][iPop]*(1.0-omega)+geqeq[iPop]*omega;

            //AntiBounce-Back conditions
//            if(iPop==0)
//                g2[xCoor][iY][iPop]=-g[xCoor2][iY][iPop];
//            else
//                g2[xCoor][iY][iPop]=-g[xCoor2][iY][iPop]+2.0*weights[iPop]*phase_temp;
        }



    }

}

void general_phase_outflow(int xCoor)
{
    int xCoor2;
    if (xCoor==0)
        xCoor2=1;
    else
        xCoor2=NX-2;

    //Update macroscopic fields
    for (int iY=1;iY<NY-1;iY++)
    {
        phase[xCoor][iY]=phase[xCoor2][iY];

		for(int iPop=0;iPop<9; iPop++)
		{

			//OutFlow conditions
			g2[xCoor][iY][iPop]=g[xCoor2][iY][iPop];

        }



    }

}

void general_density_outflow(int xCoor)
{
    int xCoor2;
    if (xCoor==0)
        xCoor2=1;
    else
        xCoor2=NX-2;

    //Update macroscopic fields
    for (int iY=1;iY<NY-1;iY++)
    {
        rho[xCoor][iY]=rho[xCoor2][iY];

		for(int iPop=0;iPop<9; iPop++)
		{

			//OutFlow conditions
			f2[xCoor][iY][iPop]=f[xCoor2][iY][iPop];

        }



    }

}



void collide_bulk()
{
    //The phase field should be calculated prior the laplacians
    for(int iX=1;iX<NX-1;iX++)
        for(int iY=1;iY<NY-1;iY++)
		{
            phase[iX][iY]=0.0;
            for(int iPop=0;iPop<9;iPop++)
   				phase[iX][iY]+=g[iX][iY][iPop];
		}

    for(int iX=1;iX<NX-1;iX++)
        for(int iY=1;iY<NY-1;iY++)
		{

			//Construction equilibrium
			rho[iX][iY]=0.0;
			ux[iX][iY]=0.0;
			uy[iX][iY]=0.0;

			for(int iPop=0;iPop<9;iPop++)
			{
				rho[iX][iY]+=f[iX][iY][iPop];
				ux[iX][iY]+=f[iX][iY][iPop]*cx[iPop];
				uy[iX][iY]+=f[iX][iY][iPop]*cy[iPop];
			}

			ux[iX][iY]=ux[iX][iY]/rho[iX][iY];
			uy[iX][iY]=uy[iX][iY]/rho[iX][iY];


			double laplace_temp=0.0;
			double gradx_temp=0.0;
			double grady_temp=0.0;
			for(int k=0;k<9;k++)
			{
				int iX2=(iX+cx[k]+NX) % NX;
				int iY2=(iY+cy[k]+NY) % NY;
				laplace_temp+=laplacestencil[k]*phase[iX2][iY2];
				gradx_temp+=gradstencilx[k]*phase[iX2][iY2];
				grady_temp+=gradstencily[k]*phase[iX2][iY2];
			}

			double phase_temp=phase[iX][iY];
			double dense_temp=rho[iX][iY];
			double ux_temp=ux[iX][iY];
			double uy_temp=uy[iX][iY];

			double sum=0.0;
			double sum_phase=0.0;
			double phase_square=phase_temp*phase_temp;
			double pressure_bulk=dense_temp/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_temp*laplace_temp;
			double chemical_pot=gammaconst*(aconst*(-phase_temp+phase_temp*phase_temp*phase_temp)-kconst*laplace_temp);

			double feqeq[9],geqeq[9];
			for (int k=1; k<9; k++)
			{
				feqeq[k]=weights[k]*(3.0*pressure_bulk+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+2.0*ux_temp*uy_temp*cx[k]*cy[k]))
				+kconst*(wxx[k]*gradx_temp*gradx_temp+wyy[k]*grady_temp*grady_temp+wxy[k]*gradx_temp*grady_temp);
				geqeq[k]=weights[k]*(3.0*chemical_pot+3.0*phase_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*phase_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+2.0*ux_temp*uy_temp*cx[k]*cy[k]));
				sum+=feqeq[k];
				sum_phase+=geqeq[k];

			}

			feqeq[0]=dense_temp-sum;
			geqeq[0]=phase_temp-sum_phase;


			for(int k=0; k < 9; k++)
			{
				f2[iX][iY][k]=f[iX][iY][k]*(1.0-omega)+omega*feqeq[k];
				g2[iX][iY][k]=g[iX][iY][k]*(1.0-omega)+omega*geqeq[k];
			}


		}

}

void update_bounce_back()
{
	//BB nodes density and velocity specification
	for(int iX=0;iX<NX;iX++)
	{
		int iXtop=(iX+1+NX)%NX;
		int iXbottom=(iX-1+NX)%NX;

		f2[iX][0][2]=f2[iX][1][4];
		f2[iX][0][5]=f2[iXtop][1][7];
		f2[iX][0][6]=f2[iXbottom][1][8];

		f2[iX][NY-1][4]=f2[iX][NY-2][2];
		f2[iX][NY-1][7]=f2[iXbottom][NY-2][5];
		f2[iX][NY-1][8]=f2[iXtop][NY-2][6];

		//BB for the scalar phase??? Ask Irina about it
		g2[iX][0][2]=g2[iX][1][4];
		g2[iX][0][5]=g2[iXtop][1][7];
		g2[iX][0][6]=g2[iXbottom][1][8];

		g2[iX][NY-1][4]=g2[iX][NY-2][2];
		g2[iX][NY-1][7]=g2[iXbottom][NY-2][5];
		g2[iX][NY-1][8]=g2[iXtop][NY-2][6];


		rho[iX][0]=1.0;
		rho[iX][NY-1]=1.0;
		phase[iX][0]=0.5;
		phase[iX][NY-1]=0.5;
		ux[iX][0]=0.0;ux[iX][NY-1]=0.0;
		uy[iX][0]=0.0;uy[iX][NY-1]=0.0;
	}

}


int main(int argc, char* argv[])
{

    matrix_init();
    init();

	//writedensity("dense_init.dat");
	//writevelocity("momentum_init.dat");


	for(int counter=0;counter<=N;counter++)
	{

        general_pressure(0,1,0,rho_inlet);
        //general_pressure(NX-1,-1,0,rho_outlet);
        //general_phase_outflow(0);
        general_phase_outflow(NX-1);
        general_density_outflow(NX-1);

        collide_bulk();

        update_bounce_back();

		//Streaming
		for(int iX=0;iX<NX;iX++)
			for(int iY=1;iY<NY-1;iY++)
				for(int iPop=0;iPop<9;iPop++)
				{
					int iX2=(iX-cx[iPop]+NX)%NX;
					int iY2=(iY-cy[iPop]+NY)%NY;
					f[iX][iY][iPop]=f2[iX2][iY2][iPop];
                    g[iX][iY][iPop]=g2[iX2][iY2][iPop];
				}

		//Writing files
		if (counter%NOUTPUT==0)
		{
			std::cout<<counter<<"\n";

			std::stringstream filewritedensity;
 			std::stringstream filewritevelocity;
 			std::stringstream filewritephase;
 			std::stringstream counterconvert;
 			counterconvert<<counter;
 			filewritedensity<<std::fixed;
			filewritevelocity<<std::fixed;
			filewritephase<<std::fixed;

			filewritedensity<<"density"<<std::string(6-counterconvert.str().size(),'0')<<counter;
			filewritevelocity<<"velocity"<<std::string(6-counterconvert.str().size(),'0')<<counter;
            filewritephase<<"phase"<<std::string(6-counterconvert.str().size(),'0')<<counter;

 			writedensity(filewritedensity.str());
			writevelocity(filewritevelocity.str());
            writephase(filewritephase.str());
		}


	}

	return 0;
}
