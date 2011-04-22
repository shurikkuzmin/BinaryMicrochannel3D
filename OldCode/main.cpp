#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkIdList.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <mpi.h>
#include "mpi_singleton.h"
#include "geometry.h"
#include "solver.h"
#include "params_list.h"

//Lattice Boltzmann initialization and parameters
const int NX=3;
const int NY=41;
const int NZ=241;
const int NPOP=19;
const int NUM=NX*NY*NZ;
const int NUMTOTAL=NUM*NPOP;
const int NUMTIME=5000;
const int NUMOUTPUT=200;
const int NUMSIGNAL=20;

//Binary-liquid initialization
const int radius=6;
const int rhol=1.0;
const int rhog=1.0;

//Parameters of the binary liquid model

int main(int argc,char* argv[])
{

	//Parallel Debugging
	#ifdef DEBUG
     	int DebugWait=1;
        while (DebugWait);
    #endif /* DEBUG */

    //Specify main communicator
	MPISingleton mpi_singleton;

    //Specify parameters
 	ParamsList params;
	params.add("aconst", 0.04);
	params.add("kconst", 0.01);
	params.add("gammaconst",1.0);
	params.add("omega",1.0);
	params.add("phase_gradient",0.5);
	params.add("phase_wall",0.0);
	params.add("rho_wall",0.5);
    params.add("tau_liq",3.0);
    params.add("tau_gas",0.7);
    params.add("tau_phi",1.0);
    params.add("force_x",0.0);
    params.add("force_y",0.0);
    params.add("force_z",0.0);
    params.add("rho_press_inlet",1.001);
    params.add("phase_press_inlet",1.0);
    params.add("rho_press_outlet",1.0);
    params.add("phase_press_outlet",1.0);


    //Specify geometry
    Geometry * geometry=new Geometry(NX,NY,NZ);
    geometry->setType(FluidNode);

    //Microchannel walls specification
//    geometry->setType(0,0,0,NX-1,0,NZ-1,SolidNode);
    //geometry->setType(0,0,0,0,NY-1,NZ-1,SolidNode);
//    geometry->setType(0,NY-1,0,NX-1,NY-1,NZ-1,SolidNode);
    //geometry->setType(NX-1,0,0,NX-1,NY-1,NZ-1,SolidNode);

    //Inlet-outlet specifications
//    geometry->setType(0,1,0,NX-1,NY-2,0,PressureNode);
//    geometry->setType(0,1,NZ-1,NX-1,NY-2,NZ-1,PressureNode);
    geometry->setType(0,0,0,NX-1,NY-1,0,PressureNode);
    geometry->setType(0,0,NZ-1,NX-1,NY-1,NZ-1,PressureNode);



    //Old Geometry
//    geometry->setType(0,0,0,NX-1,NY-1,0,SolidNode);
//    geometry->setType(0,0,NZ-1,NX-1,NY-1,NZ-1,SolidNode);
//    geometry->setType(0,0,0,NX-1,0,NZ-1,SolidNode);
//    geometry->setType(0,NY-1,0,NX-1,NY-1,NZ-1,SolidNode);
    //Specify solver



    Solver solver(geometry,params);

	//Density and phase fields initialization
	double rho_temp;
	double phase_temp;
	double u_temp[3];
	for (int counter=0; counter<NUM; counter++)
	{
		int iZ=counter/(NX*NY);
		int iY=(counter%(NX*NY))/NX;
		int iX=(counter%(NX*NY))%NX;

		//Initialization for the droplets coalescence
//		if(((iX-floor(NX/4))*(iX-floor(NX/4))+(iY-floor(NY/2))*(iY-floor(NY/2))+(iZ-floor(NZ/2))*(iZ-floor(NZ/2))<=radius*radius) ||
//		   ((iX-floor(3*NX/4))*(iX-floor(3*NX/4))+(iY-floor(NY/2))*(iY-floor(NY/2))+(iZ-floor(NZ/2))*(iZ-floor(NZ/2))<=radius*radius))
//		{
//			rho_temp=rhol;
//			phase_temp=1.0;
//            if(iX<NX/2)
//                u_temp[0]=0.02;
//            else
//                u_temp[0]=-0.02;
//                u_temp[1]=0.0;
//                u_temp[2]=0.0;
//
//        }
//		else
//		{
//			rho_temp=rhog;
//			phase_temp=-1.0;
//            u_temp[0]=0.0;
//            u_temp[1]=0.0;
//            u_temp[2]=0.0;
//
//		}


        //Initialization for the droplet on the surface
//		if ((iX-floor(NX/2))*(iX-floor(NX/2))+(iY-floor(NY/2))*(iY-floor(NY/2))+(iZ-floor(radius))*(iZ-floor(radius))<=radius*radius)
//		{
//			rho_temp=rhol;
//			phase_temp=1.0;
//            u_temp[0]=0.0;
//            u_temp[1]=0.0;
//            u_temp[2]=0.0;
//
//        }
//		else
//        {
//			rho_temp=rhog;
//			phase_temp=-1.0;
//            u_temp[0]=0.0;
//            u_temp[1]=0.0;
//            u_temp[2]=0.0;
//
//		}
//

        //Initialization of the part of the channel

//		if ((iX>=NX/3)&&(iX<=2*NX/3)&&(iY>=3)&&(iY<NY-3)&&(iZ>=3)&&(iZ<NZ-3))
//		{
//			rho_temp=rhol;
//			phase_temp=1.0;
//            u_temp[0]=0.0;
//            u_temp[1]=0.0;
//            u_temp[2]=0.0;
//
//        }
//		else
//        {
//			rho_temp=rhog;
//			phase_temp=-1.0;
//            u_temp[0]=0.0;
//            u_temp[1]=0.0;
//            u_temp[2]=0.0;
//
//		}

		if ((iZ>=NZ/6)&&(iZ<=5*NZ/6)&&(iY>=4)&&(iY<NY-4))
		{
			rho_temp=rhog;
			phase_temp=-1.0;
            u_temp[0]=0.0;
            u_temp[1]=0.0;
            u_temp[2]=0.0;

        }
		else
        {
			rho_temp=rhol;
			phase_temp=1.0;
            u_temp[0]=0.0;
            u_temp[1]=0.0;
            u_temp[2]=0.0;

		}





		solver.putDensity(iX,iY,iZ,iX,iY,iZ,rho_temp);
		solver.putPhase(iX,iY,iZ,iX,iY,iZ,phase_temp);
		solver.putVelocity(iX,iY,iZ,iX,iY,iZ,u_temp);
    }

    //Initialization of the populations
    solver.init();

	//Main iteration
	for (int time_counter=0; time_counter<NUMTIME; time_counter++)
	{
		//Collision procedure with the reconstruction of the equilibrium populations
        solver.collide_stream();

		if (time_counter%NUMSIGNAL==0)
		{
            cout<<"Time is "<<time_counter<<"\n";
        }

		//Output files
		if (time_counter%NUMOUTPUT==0)
		{
            std::stringstream file_name;
   			std::stringstream time_string;
 			time_string<<time_counter;

			file_name<<"phase"<<std::string(5-time_string.str().size(),'0')<<time_counter;
            //solver.writeWholeDensityPhaseVelocity(file_name.str());
            solver.writeTextYZVelocity(file_name.str());
            cout<<"The iteration step is "<<time_counter<<"\n";
		}
	}

    //Destroy objects
    delete geometry;

	return 0;
}
