#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

const int NY=8;
//const int NX=(NY-2)*16; 
const int NX=11;
const int N=100000;	

double f[NX][NY][9], f2[NX][NY][9];
double rho[NX][NY],ux[NX][NY],uy[NX][NY];

double rho1=1.003;
double rho2=1.0;

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



int main(int argc, char* argv[])
{

	double omega=1.0;
	//double forcex=0.0001/1.0*(2.0-omega)/omega;
	//double forcey=0.000;

	double omegaginzburg=8.0*(2.0-omega)/(8.0-omega);
	double omegamat[]={1.0,1.0,1.0,omega,omega,omega,1.0,omegaginzburg,omegaginzburg};


	
	double weights[]={4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};
	int cx[]={0,1,0,-1,0,1,-1,-1,1};
	int cy[]={0,0,1,0,-1,1,1,-1,-1};
	double M[9][9];
	
	double g[]={1.0,-2.0,-2.0,-2.0,-2.0,4.0,4.0,4.0,4.0};
	for (int iCoor=0;iCoor<9;iCoor++)
	{
		M[0][iCoor]=1.0;
		M[1][iCoor]=cx[iCoor]*sqrt(3.0);
		M[2][iCoor]=cy[iCoor]*sqrt(3.0);
		M[3][iCoor]=(cx[iCoor]*cx[iCoor]-1.0/3.0)*3.0/sqrt(2.0);
		M[4][iCoor]=cx[iCoor]*cy[iCoor]*3.0;
		M[5][iCoor]=(cy[iCoor]*cy[iCoor]-1.0/3.0)*3.0/sqrt(2.0);
		M[6][iCoor]=g[iCoor]/2.0;
		M[7][iCoor]=g[iCoor]*cx[iCoor]*sqrt(1.5)/2.0;
		M[8][iCoor]=g[iCoor]*cy[iCoor]*sqrt(1.5)/2.0;
	}
	
	
	//Initialization of initial conditions
	for(int iX=1;iX<NX-1;iX++)
		for(int iY=0;iY<NY;iY++)
		{
			rho[iX][iY]=1.0;
			ux[iX][iY]=0.0;
			uy[iX][iY]=0.0;
			
			double fluxx=3.0*ux[iX][iY];
			double fluxy=3.0*uy[iX][iY];
			
			double qxx=4.5*ux[iX][iY]*ux[iX][iY];
			double qxy=9.0*ux[iX][iY]*uy[iX][iY];
			double qyy=4.5*uy[iX][iY]*uy[iX][iY];
	
			for(int iPop=0;iPop<9;iPop++)
			{
				f[iX][iY][iPop]=weights[iPop]*(rho[iX][iY]+fluxx*cx[iPop]+fluxy*cy[iPop]+
								qxx*(cx[iPop]*cx[iPop]-1.0/3.0)+qxy*cx[iPop]*cy[iPop]+qyy*(cy[iPop]*cy[iPop]-1.0/3.0));	
			}
		}
	for(int iY=0;iY<NY;iY++)
	{
		rho[0][iY]=rho1;
		ux[0][iY]=0.0;
		uy[0][iY]=0.0;

		rho[NX-1][iY]=rho2;
		ux[NX-1][iY]=0.0;
		uy[NX-1][iY]=0.0;
		
	
		for(int iPop=0;iPop<9;iPop++)
		{
			f[0][iY][iPop]=weights[iPop]*rho[0][iY];	
			f[NX-1][iY][iPop]=weights[iPop]*rho[NX-1][iY];	

		}
	}


	for(int counter=0;counter<=N;counter++)
	{
		for(int iX=0;iX<NX;iX++)
			for(int iY=0;iY<NY;iY++)
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
				
				double dense=rho[iX][iY];

				//ux[iX][iY]/=dense;
				//uy[iX][iY]/=dense;
				
				//ux[iX][iY]+=forcex/(2.0*dense);
				//ux[iX][iY]+=forcey/(2.0*dense);
				
				double uxeq=ux[iX][iY];
				double uyeq=uy[iX][iY];
	
				//Construction of the equilibrium moments
				double eq[9];
				eq[0]=dense;
				eq[1]=sqrt(3.0)*dense*uxeq;
				eq[2]=sqrt(3.0)*dense*uyeq;
				eq[3]=3.0/sqrt(2.0)*dense*uxeq*uxeq;
				eq[4]=3.0*dense*uxeq*uyeq;
				eq[5]=3.0/sqrt(2.0)*dense*uyeq*uyeq;
				eq[6]=4.5*dense*uxeq*uxeq*uyeq*uyeq;
				eq[7]=3.0*sqrt(1.5)*dense*uxeq*uyeq*uyeq;
				eq[8]=3.0*sqrt(1.5)*dense*uxeq*uxeq*uyeq;
	
				//double feq;
				//double feqforce;
				
				double feqeq[9];	
				for(int iPop=0;iPop<9; iPop++)
				{	
// 					feqeq[iPop]=weights[iPop]*dense*(1.0+3.0*cx[iPop]*uxeq+3.0*cy[iPop]*uyeq+
// 					4.5*(cx[iPop]*cx[iPop]-1.0/3.0)*uxeq*uxeq+9.0*cx[iPop]*cy[iPop]*uxeq*uyeq+4.5*(cy[iPop]*cy[iPop]-1.0/3.0)*uyeq*uyeq+
// 					81.0/4.0*(cx[iPop]*cx[iPop]*cy[iPop]*cy[iPop]-1.0/3.0*cx[iPop]*cx[iPop]-1.0/3.0*cy[iPop]*cy[iPop]+1.0/9.0)*uxeq*uxeq*uyeq*uyeq+
// 					27.0/2.0*(cx[iPop]*cx[iPop]*cy[iPop]-1.0/3.0*cy[iPop])*uxeq*uxeq*uyeq+
// 					27.0/2.0*(cy[iPop]*cy[iPop]*cx[iPop]-1.0/3.0*cx[iPop])*uyeq*uyeq*uxeq);
				
					feqeq[iPop]=weights[iPop]*(dense+3.0*cx[iPop]*uxeq+3.0*cy[iPop]*uyeq+
						4.5*(cx[iPop]*cx[iPop]-1.0/3.0)*uxeq*uxeq+9.0*cx[iPop]*cy[iPop]*uxeq*uyeq+4.5*(cy[iPop]*cy[iPop]-1.0/3.0)*uyeq*uyeq);
				}
				
				double add[9];
				double addit;			
				for(int iPop=0;iPop < 9; iPop++)
				{
					add[iPop]=0.0;
					for(int k=0; k < 9; k++)
						add[iPop]=add[iPop]+M[iPop][k]*(-f[iX][iY][k]); //+feqeq[k]); //+node.popeq[k]);
					add[iPop]+=eq[iPop];
				}
					
				for(int k=0; k < 9; k++)
				{
 					//feqforce=(1.0-0.5*omega)*weights[k]*(forcex*(3.0*(cx[k]-uxeq)+9.0*cx[k]*(cx[k]*uxeq+cy[k]*uyeq))+
					//		forcey*(3.0*(cy[k]-uyeq)+9.0*cy[k]*(cx[k]*uxeq+cy[k]*uyeq)));
					
 					//feqforce=weights[k]*(3.0*(1.0-0.5*omegamat[1])*forcex*cx[k]+3.0*(1.0-0.5*omegamat[2])*forcey*cy[k]+
 					/*		9.0*((1.0-0.5*omegamat[3])*(cx[k]*cx[k]-1.0/3.0)*forcex*uxeq+
 							(1.0-0.5*omegamat[4])*cx[k]*cy[k]*(forcex*uyeq+forcey*uxeq)+
 							(1.0-0.5*omegamat[5])*(cy[k]*cy[k]-1.0/3.0)*forcey*uyeq));	
					*/	
					
					addit=0.0;
					for(int m=0; m < 9; m++)
						addit=addit+omegamat[m]*M[m][k]*add[m];
					//f2[iX][iY][k]=f[iX][iY][k]+weights[k]*addit; //+feqforce;
					f2[iX][iY][k]=f[iX][iY][k]*(1.0-omega)+omega*feqeq[k]; //+feqforce;
				}
					
				
			}
		
		//Boundary conditions and BB nodes density and velocity specification
/*		for(int iX=0;iX<NX;iX++)
		{
			int iXtop=(iX+1+NX)%NX;
			int iXbottom=(iX-1+NX)%NX;
			
			f2[iX][0][2]=f2[iX][1][4];
			f2[iX][0][5]=f2[iXtop][1][7];
			f2[iX][0][6]=f2[iXbottom][1][8];
			
			f2[iX][NY-1][4]=f2[iX][NY-2][2];
			f2[iX][NY-1][7]=f2[iXbottom][NY-2][5];
			f2[iX][NY-1][8]=f2[iXtop][NY-2][6];
			
			rho[iX][0]=1.0;
			rho[iX][NY-1]=1.0;
			ux[iX][0]=0.0;ux[iX][NY-1]=0.0;
			uy[iX][0]=0.0;uy[iX][NY-1]=0.0;
		}		*/
	
		
		//Streaming
		for(int iX=0;iX<NX;iX++)
			for(int iY=0;iY<NY;iY++)
				for(int iPop=0;iPop<9;iPop++)
				{
					int iX2=(iX-cx[iPop]+NX)%NX;
					int iY2=(iY-cy[iPop]+NY)%NY;
					f[iX][iY][iPop]=f2[iX2][iY2][iPop];
				}
		//Pressure boundary conditions
		for (int iY=1;iY<NY-1;iY++)	
		{
			ux[0][iY]=1.-(f[0][iY][0]+f[0][iY][2]+f[0][iY][4]+2.0*(f[0][iY][3]+f[0][iY][6]+f[0][iY][7]))/rho1;
			f[0][iY][1] = f[0][iY][3]+2.0*rho1*ux[0][iY]/3.0;
			f[0][iY][5] = f[0][iY][7]+0.5*(f[0][iY][4]-f[0][iY][2])+rho1*ux[0][iY]/6.;
			f[0][iY][8] = f[0][iY][6]-0.5*(f[0][iY][4]-f[0][iY][2])+rho1*ux[0][iY]/6.;
			// East Boundary :
			ux[NX-1][iY]=-1.+(f[NX-1][iY][0]+f[NX-1][iY][2]+f[NX-1][iY][4]+2.*(f[NX-1][iY][1]+f[NX-1][iY][5]+f[NX-1][iY][8]))/rho2;
			f[NX-1][iY][3] = f[NX-1][iY][1]-2.0*rho2*ux[NX-1][iY]/3.;
			f[NX-1][iY][7] = f[NX-1][iY][5]-0.5*(f[NX-1][iY][4]-f[NX-1][iY][2])-rho2*ux[NX-1][iY]/6.;
			f[NX-1][iY][6] = f[NX-1][iY][8]+0.5*(f[NX-1][iY][4]-f[NX-1][iY][2])-rho2*ux[NX-1][iY]/6.;
		}		
		
		//Wall Zou-He conditions with zero
		for (int iX=1;iX<NX-1;iX++)	
		{
			f[iX][0][2] = f[iX][0][4];
			f[iX][0][5] = f[iX][0][7];
			f[iX][0][6] = f[iX][0][8];
			// East Boundary :
			f[iX][NY-1][4] = f[iX][NY-1][2];
			f[iX][NY-1][7] = f[iX][NY-1][5];
			f[iX][NY-1][8] = f[iX][NY-1][6];
		}		
		//Corners
		f[0][0][1]=f[0][0][3];
		f[0][0][2]=f[0][0][4];
		f[0][0][5]=f[0][0][7];
		f[0][0][6]=0.5*(rho1-(f[0][0][0]+f[0][0][1]+f[0][0][2]+f[0][0][3]+f[0][0][4]+f[0][0][5]+f[0][0][7]));
		f[0][0][8]=f[0][0][6];

		f[0][NY-1][1]=f[0][NY-1][3];
		f[0][NY-1][4]=f[0][NY-1][2];
		f[0][NY-1][8]=f[0][NY-1][6];
		f[0][NY-1][5]=0.5*(rho1-(f[0][NY-1][0]+f[0][NY-1][1]+f[0][NY-1][2]+f[0][NY-1][3]+f[0][NY-1][4]+f[0][NY-1][6]+f[0][NY-1][8]));
		f[0][NY-1][7]=f[0][NY-1][5];

		f[NX-1][0][3]=f[NX-1][0][1];
		f[NX-1][0][2]=f[NX-1][0][4];
		f[NX-1][0][6]=f[NX-1][0][8];
		f[NX-1][0][5]=0.5*(rho2-(f[NX-1][0][0]+f[NX-1][0][1]+f[NX-1][0][2]+f[NX-1][0][3]+f[NX-1][0][4]+f[NX-1][0][6]+f[NX-1][0][8]));
		f[NX-1][0][7]=f[NX-1][0][5];

		f[NX-1][NY-1][3]=f[NX-1][NY-1][1];
		f[NX-1][NY-1][4]=f[NX-1][NY-1][2];
		f[NX-1][NY-1][7]=f[NX-1][NY-1][5];
		f[NX-1][NY-1][6]=0.5*(rho2-(f[NX-1][NY-1][0]+f[NX-1][NY-1][1]+f[NX-1][NY-1][2]+f[NX-1][NY-1][3]+f[NX-1][NY-1][4]+f[NX-1][NY-1][5]+f[NX-1][NY-1][7]));
		f[NX-1][NY-1][8]=f[NX-1][NY-1][6];

		//Writing files
		if (counter%1000==0)
		{
			std::cout<<counter<<"\n";
			
			std::stringstream filewritedensity;
 			std::stringstream filewritevelocity;
 			std::stringstream counterconvert;
 			counterconvert<<counter;
 			filewritedensity<<std::fixed;
			filewritevelocity<<std::fixed;
			
			filewritedensity<<"proba"<<std::string(6-counterconvert.str().size(),'0')<<counter;
			filewritevelocity<<"phase"<<std::string(6-counterconvert.str().size(),'0')<<counter;
 			
 			writedensity(filewritedensity.str());
			writevelocity(filewritevelocity.str());	
		}
		
	
	}

	return 0;
}