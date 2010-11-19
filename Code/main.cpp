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
const int NX=30;
const int NY=30;
const int NZ=451;
const int NPOP=19;
const int NUM=NX*NY*NZ;
const int NUMTOTAL=NUM*NPOP;
const int NUMTIME=10000;
const int NUMOUTPUT=1000;
const int NUMSIGNAL=20;

//Binary-liquid initialization
const int width=4;
const int radius=6;
const double rhol=1.0;
const double rhog=1.0;

void macro_init(Solver & solver)
{
    //Density and phase field initialization
	double rho_temp;
	double phase_temp;
	double u_temp[3];
	for (int counter=0; counter<NUM; counter++)
	{
		int iZ=counter/(NX*NY);
		int iY=(counter%(NX*NY))/NX;
		int iX=(counter%(NX*NY))%NX;

        //Initialization of the part of the channel

		if ((iZ>=(NZ-1)/3)&&(iZ<=2*(NZ-1)/3)&&(iX>=width)&&(iX<=NX-width-1)&&(iY>=width)&&(iY<=NY-width-1))
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

}

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
	params.add("phase_gradient",0.0);
	params.add("phase_wall",0.0);
	params.add("rho_wall",0.5);
    params.add("tau_liq",2.7);
    params.add("tau_gas",0.7);
    params.add("force_x",0.0);
    params.add("force_y",0.0);
    params.add("force_z",0.0001);

    //Useless parameters to be deleted after
    params.add("rho_press_inlet",1.0);
    params.add("rho_press_outlet",1.0);
    params.add("phase_press_inlet",1.0);
    params.add("phase_press_outlet",1.0);

    //Specify geometry
    Geometry * geometry=new Geometry(NX,NY,NZ);
    geometry->setType(FluidNode);

    //Microchannel walls specification
    geometry->setType(0,0,0,NX-1,0,NZ-1,SolidNode);
    geometry->setType(0,0,0,0,NY-1,NZ-1,SolidNode);
    geometry->setType(0,NY-1,0,NX-1,NY-1,NZ-1,SolidNode);
    geometry->setType(NX-1,0,0,NX-1,NY-1,NZ-1,SolidNode);

    Solver solver(geometry,params);

    //Solver initialization from the file
    solver.load_file("equili");

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
            if (solver.checkNAN())
            {
                cout<<"Phase fields contain NaN values\n";
                MPI_Abort(MPI_COMM_WORLD,-1);
            }

		}
		//Output files
		if (time_counter%NUMOUTPUT==0)
		{
            std::stringstream file_name;
   			std::stringstream time_string;
 			time_string<<time_counter;

			file_name<<"phase"<<std::string(5-time_string.str().size(),'0')<<time_counter;
            solver.writeWholeDensityPhaseVelocity(file_name.str());
            cout<<"Output is done on the step "<<time_counter<<"\n";
		}
	}

    //Destroy objects
    delete geometry;
	return 0;
}
