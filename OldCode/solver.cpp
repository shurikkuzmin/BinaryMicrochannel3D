#include <vtkXMLStructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkIdList.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include "solver.h"
#include "dynamics_BGK.h"
#include "dynamics_special.h"
#include "dynamics_simple_BGK.h"
#include "dynamics_simple_special.h"
#include "dynamics_none.h"
#include "dynamics_zouhe_pressure.h"
#include "dynamics_zouhe_simple_pressure.h"
#include <mpi.h>
#include <iostream>
#include <sstream>
#include <string>

const char Solver::cx[19]={0,1,-1,0, 0,0, 0,1,-1, 1,-1,0, 0, 0, 0,1,-1, 1,-1};
const char Solver::cy[19]={0,0, 0,1,-1,0, 0,1, 1,-1,-1,1,-1, 1,-1,0, 0, 0, 0};
const char Solver::cz[19]={0,0, 0,0, 0,1,-1,0, 0, 0, 0,1, 1,-1,-1,1, 1,-1,-1};

Solver::Solver(Geometry * _geom,ParamsList _params_list):
	geom(_geom),NX(geom->getNX()),
	NY(geom->getNY()),NZ(geom->getNZ()),
	rank(mpi_wrapper().get_rank()),
	size(mpi_wrapper().get_size()),
	params_list(_params_list),
	NPOP(19)
{
	//Initialization of ranges
	int zstep=NZ/size;
	int zbegin=rank*zstep;
	int zend=rank*zstep+zstep-1;
	if(rank==size-1)
		zend=NZ-1;
	int NUMLOCAL=NX*NY*(zend-zbegin+1);

	//Creation the necessary objects
	lattice = new Lattice(zbegin,zend,NX,NY,NZ,NUMLOCAL);
	Dynamics * dynamics_bgk;
	#ifdef HYDRO
        dynamics_bgk = destroyer_list.take_ownership(new DynamicsSimpleBGK(this,params_list));
	#else
        dynamics_bgk = destroyer_list.take_ownership(new DynamicsBGK(this,params_list));
    #endif
    Dynamics * dynamics_special;
    Dynamics * dynamics_none = destroyer_list.take_ownership(new DynamicsNONE(this,params_list));

    //At this moment we just assume that pressure boundaries are on Z location


    Dynamics * dynamics_pressure_inlet;
    Dynamics * dynamics_pressure_outlet;

	params_list.add("rho_press",params_list("rho_press_inlet").value<double>());
	params_list.add("phase_press",params_list("phase_press_inlet").value<double>());
    params_list.add("normx",0);
    params_list.add("normy",0);
    params_list.add("normz",1);

    #ifdef HYDRO
        dynamics_pressure_inlet=destroyer_list.take_ownership(new DynamicsZouHeSimplePressure(this,params_list));
    #else
        dynamics_pressure_inlet=destroyer_list.take_ownership(new DynamicsZouHePressure(this,params_list));
    #endif
    params_list.change("rho_press",params_list("rho_press_outlet").value<double>());
    params_list.change("phase_press",params_list("phase_press_outlet").value<double>());
    params_list.change("normz",-1);
    #ifdef HYDRO
        dynamics_pressure_outlet=destroyer_list.take_ownership(new DynamicsZouHeSimplePressure(this,params_list));
    #else
        dynamics_pressure_outlet=destroyer_list.take_ownership(new DynamicsZouHePressure(this,params_list));
    #endif
	//Allocation of memory
	if (rank==0)
	{
		density=new double[NX*NY*NZ];
		phase = new double[NX*NY*NZ];
        velocity=new double[3*NX*NY*NZ];
	}


	for (int iZ=0;iZ<NZ;iZ++)
		for(int iY=0;iY<NY;iY++)
			for(int iX=0;iX<NX;iX++)
			{
				if ((iZ>=zbegin)&&(iZ<=zend))
				{
                    lattice->geometry[(iZ-zbegin)*NX*NY+iY*NX+iX]=geom->getType(iX,iY,iZ);
                    if (geom->getType(iX,iY,iZ)==FluidNode)
                    {
                        bool flag=false;

                        for(int k=0;k<NPOP;k++)
                        {
                            int iX2=(iX+cx[k]+NX)%NX;
                            int iY2=(iY+cy[k]+NY)%NY;
                            int iZ2=(iZ+cz[k]+NZ)%NZ;
                            if (geom->getType(iX2,iY2,iZ2)==SolidNode)
                                flag=true;
                        }
                        if (flag==true)
                        {
                            iterX=iX;
                            iterY=iY;
                            iterZ=iZ;
                            #ifdef HYDRO
                                dynamics_special=destroyer_list.take_ownership(new DynamicsSimpleSpecial(this,params_list));
                            #else
                                dynamics_special=destroyer_list.take_ownership(new DynamicsSpecial(this,params_list));
                            #endif
                            lattice->dynamics_list[(iZ-zbegin)*NX*NY+iY*NX+iX]=dynamics_special;
                        }
                        else
                            lattice->dynamics_list[(iZ-zbegin)*NX*NY+iY*NX+iX]=dynamics_bgk;
                    }
       				else if (geom->getType(iX,iY,iZ)==PressureNode)
       				{
                        //We do it simple without derivation of normals right now
                        if (iZ==0)
                        {
                            lattice->dynamics_list[(iZ-zbegin)*NX*NY+iY*NX+iX]=dynamics_pressure_inlet;
                        }
                        else if (iZ==NZ-1)
                        {
                            lattice->dynamics_list[(iZ-zbegin)*NX*NY+iY*NX+iX]=dynamics_pressure_outlet;
                        }

                    }
       				else
                        lattice->dynamics_list[(iZ-zbegin)*NX*NY+iY*NX+iX]=dynamics_none;

				}
			}

    //Preparation of the geometry types above and at the bottom
	int zbottom=(zbegin-1+NZ)%NZ;
	int ztop=(zend+1+NZ)%NZ;
	for(int iY=0;iY<NY;iY++)
		for(int iX=0;iX<NX;iX++)
		{
            lattice->geom_top[iY*NX+iX]=geom->getType(iX,iY,ztop);
            lattice->geom_bottom[iY*NX+iX]=geom->getType(iX,iY,zbottom);
		}

}


void Solver::putDensity(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,double dense)
{
    int zbegin=lattice->getZbegin();
    int zend=lattice->getZend();
    if (iZbegin>zend)
        return;
    if (iZend<zbegin)
        return;
    if (iZend>zend) iZend=zend;
    if (iZbegin<zbegin) iZbegin=zbegin;

    //cout<<"zend="<<zend<<" zbegin="<<zbegin;
    lattice->putDensity(iXbegin,iYbegin,iZbegin-zbegin,iXend,iYend,iZend-zbegin,dense);
}
void Solver::putDensity(double dense)
{
    int zbegin=lattice->getZbegin();
    int zend=lattice->getZend();
    lattice->putDensity(0,0,0,NX-1,NY-1,zend-zbegin,dense);
}

void Solver::putPhase(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,double phase)
{
    int zbegin=lattice->getZbegin();
    int zend=lattice->getZend();
    if (iZbegin>zend)
        return;
    if (iZend<zbegin)
        return;
    if (iZend>zend) iZend=zend;
    if (iZbegin<zbegin) iZbegin=zbegin;

    lattice->putPhase(iXbegin,iYbegin,iZbegin-zbegin,iXend,iYend,iZend-zbegin,phase);
}
void Solver::putPhase(double phase)
{
    int zbegin=lattice->getZbegin();
    int zend=lattice->getZend();
    lattice->putPhase(0,0,0,NX-1,NY-1,zend-zbegin,phase);
}

void Solver::putVelocity(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,double * velocity)
{
    int zbegin=lattice->getZbegin();
    int zend=lattice->getZend();
    if (iZbegin>zend)
        return;
    if (iZend<zbegin)
        return;
    if (iZend>zend) iZend=zend;
    if (iZbegin<zbegin) iZbegin=zbegin;

    lattice->putVelocity(iXbegin,iYbegin,iZbegin-zbegin,iXend,iYend,iZend-zbegin,velocity);
}
void Solver::putVelocity(double * velocity)
{
    int zbegin=lattice->getZbegin();
    int zend=lattice->getZend();
    lattice->putVelocity(0,0,0,NX-1,NY-1,zend-zbegin,velocity);
}


void Solver::getWholeDensity()
{
	if (size==1)
	{
		for(int iCount=0;iCount<NX*NY*NZ;iCount++)
			density[iCount]=lattice->rho[iCount];
		return;
	}
    else
	{
		if (rank==0)
		{
       		for(int iCount=0;iCount<lattice->getNUMLOCAL();iCount++)
                density[iCount]=lattice->rho[iCount];

			int offset=lattice->getNUMLOCAL();
			for (int nproc=1;nproc<size;nproc++)
			{

				MPI_Status status;
				int sizes[2];

				MPI_Recv(sizes,2,MPI_INT,nproc,2,MPI_COMM_WORLD,&status);

				double * temp_array=new double[(sizes[1]-sizes[0]+1)*NX*NY];
				MPI_Recv(temp_array,(sizes[1]-sizes[0]+1)*NX*NY,MPI_DOUBLE,nproc,3,MPI_COMM_WORLD,&status);

 				for(int iCount=0;iCount<(sizes[1]-sizes[0]+1)*NX*NY;iCount++)
					density[iCount+offset]=temp_array[iCount];

				offset=offset+(sizes[1]-sizes[0]+1)*NX*NY;

				delete[] temp_array;
			}
		}
		else
		{
			int sizes[2]={lattice->getZbegin(),lattice->getZend()};
			MPI_Send(sizes,2,MPI_INT,0,2,MPI_COMM_WORLD);
			MPI_Send(lattice->rho,(lattice->getZend()-lattice->getZbegin()+1)*NX*NY,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);

	}
}

void Solver::getWholePhase()
{
	if (size==1)
	{
		for(int iCount=0;iCount<NX*NY*NZ;iCount++)
			phase[iCount]=lattice->phase[iCount];
		return;
	}
    else
	{
		if (rank==0)
		{
       		for(int iCount=0;iCount<lattice->getNUMLOCAL();iCount++)
                phase[iCount]=lattice->phase[iCount];

			int offset=lattice->getNUMLOCAL();
			for (int nproc=1;nproc<size;nproc++)
			{

				MPI_Status status;
				int sizes[2];

				MPI_Recv(sizes,2,MPI_INT,nproc,2,MPI_COMM_WORLD,&status);

				double * temp_array=new double[(sizes[1]-sizes[0]+1)*NX*NY];
				MPI_Recv(temp_array,(sizes[1]-sizes[0]+1)*NX*NY,MPI_DOUBLE,nproc,3,MPI_COMM_WORLD,&status);

 				for(int iCount=0;iCount<(sizes[1]-sizes[0]+1)*NX*NY;iCount++)
					phase[iCount+offset]=temp_array[iCount];

				offset=offset+(sizes[1]-sizes[0]+1)*NX*NY;

				delete[] temp_array;
			}
		}
		else
		{
			int sizes[2]={lattice->getZbegin(),lattice->getZend()};
			MPI_Send(sizes,2,MPI_INT,0,2,MPI_COMM_WORLD);
			MPI_Send(lattice->phase,(lattice->getZend()-lattice->getZbegin()+1)*NX*NY,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);

	}
}

void Solver::getWholeVelocity()
{
	if (size==1)
	{
		for(int iCount=0;iCount<NX*NY*NZ;iCount++)
        {
			velocity[3*iCount]=lattice->ux[iCount];
            velocity[3*iCount+1]=lattice->uy[iCount];
            velocity[3*iCount+2]=lattice->uz[iCount];
        }
		return;
	}
    else
	{
		if (rank==0)
		{
       		for(int iCount=0;iCount<lattice->getNUMLOCAL();iCount++)
            {
                velocity[3*iCount]=lattice->ux[iCount];
                velocity[3*iCount+1]=lattice->uy[iCount];
                velocity[3*iCount+2]=lattice->uz[iCount];
            }

			int offset=3*lattice->getNUMLOCAL();
			for (int nproc=1;nproc<size;nproc++)
			{

				MPI_Status status;
				int sizes[2];

				MPI_Recv(sizes,2,MPI_INT,nproc,2,MPI_COMM_WORLD,&status);

				double * temp_array_x=new double[(sizes[1]-sizes[0]+1)*NX*NY];
				double * temp_array_y=new double[(sizes[1]-sizes[0]+1)*NX*NY];
				double * temp_array_z=new double[(sizes[1]-sizes[0]+1)*NX*NY];
				MPI_Recv(temp_array_x,(sizes[1]-sizes[0]+1)*NX*NY,MPI_DOUBLE,nproc,3,MPI_COMM_WORLD,&status);
				MPI_Recv(temp_array_y,(sizes[1]-sizes[0]+1)*NX*NY,MPI_DOUBLE,nproc,3,MPI_COMM_WORLD,&status);
				MPI_Recv(temp_array_z,(sizes[1]-sizes[0]+1)*NX*NY,MPI_DOUBLE,nproc,3,MPI_COMM_WORLD,&status);


 				for(int iCount=0;iCount<(sizes[1]-sizes[0]+1)*NX*NY;iCount++)
				{
					velocity[3*iCount+offset]=temp_array_x[iCount];
					velocity[3*iCount+1+offset]=temp_array_y[iCount];
					velocity[3*iCount+2+offset]=temp_array_z[iCount];
				}
				offset=offset+3*(sizes[1]-sizes[0]+1)*NX*NY;

				delete[] temp_array_x;
				delete[] temp_array_y;
				delete[] temp_array_z;
			}
		}
		else
		{
			int sizes[2]={lattice->getZbegin(),lattice->getZend()};
			MPI_Send(sizes,2,MPI_INT,0,2,MPI_COMM_WORLD);
			MPI_Send(lattice->ux,(lattice->getZend()-lattice->getZbegin()+1)*NX*NY,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
			MPI_Send(lattice->uy,(lattice->getZend()-lattice->getZbegin()+1)*NX*NY,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
			MPI_Send(lattice->uz,(lattice->getZend()-lattice->getZbegin()+1)*NX*NY,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);

	}
}




void Solver::writeWholeDensity(std::string name)
{
	getWholeDensity();

	if(rank==0)
	{
        name=name+".vts";

        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                for(int counterX=0;counterX<NX;counterX++)
                    points->InsertNextPoint(counterX,counterY,counterZ);


        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfComponents(1);
        data->SetName("Density");



        for(vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
        {
            double * temp_data=&density[i];
            data->InsertNextTupleValue(temp_data);
        }

        vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
        structuredGrid->SetDimensions(NX,NY,NZ);
        structuredGrid->SetPoints(points);
        structuredGrid->GetPointData()->AddArray(data);

        vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
        writer->SetFileName(name.c_str());
        //writer.SetPoints();
        writer->SetInput(structuredGrid);
        writer->Write();
	}
}

void Solver::writeWholePhase(std::string name)
{
    Solver::getWholePhase();
    if (rank==0)
    {
        name=name+".vts";

        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                for(int counterX=0;counterX<NX;counterX++)
                    points->InsertNextPoint(counterX,counterY,counterZ);


        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfComponents(1);
        data->SetName("Phase");



        for(vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
        {
            double * temp_data=&phase[i];
            data->InsertNextTupleValue(temp_data);
        }

        vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
        structuredGrid->SetDimensions(NX,NY,NZ);
        structuredGrid->SetPoints(points);
        structuredGrid->GetPointData()->AddArray(data);

        vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
        writer->SetFileName(name.c_str());
        //writer.SetPoints();
        writer->SetInput(structuredGrid);
        writer->Write();
    }
}

void Solver::writeWholeVelocity(std::string name)
{
    Solver::getWholeVelocity();
    if (rank==0)
    {
        name=name+".vts";

        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                for(int counterX=0;counterX<NX;counterX++)
                    points->InsertNextPoint(counterX,counterY,counterZ);


        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfComponents(3);
        data->SetName("Velocity");



        for(vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
        {
            double * temp_data=&velocity[3*i];
            data->InsertNextTupleValue(temp_data);
        }

        vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
        structuredGrid->SetDimensions(NX,NY,NZ);
        structuredGrid->SetPoints(points);
        structuredGrid->GetPointData()->AddArray(data);

        vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
        writer->SetFileName(name.c_str());
        writer->SetInput(structuredGrid);
        writer->Write();
    }
}

void Solver::writeWholeDensityPhaseVelocity(std::string name)
{

    Solver::getWholeDensity();
    Solver::getWholePhase();
    Solver::getWholeVelocity();
    if (rank==0)
    {
        name=name+".vts";

        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                for(int counterX=0;counterX<NX;counterX++)
                    points->InsertNextPoint(counterX,counterY,counterZ);

        vtkSmartPointer<vtkDoubleArray> data_phase = vtkSmartPointer<vtkDoubleArray>::New();
        vtkSmartPointer<vtkDoubleArray> data_density = vtkSmartPointer<vtkDoubleArray>::New();
        vtkSmartPointer<vtkDoubleArray> data_velocity = vtkSmartPointer<vtkDoubleArray>::New();

        data_phase->SetNumberOfComponents(1);
        data_phase->SetName("Phase");
        data_density->SetNumberOfComponents(1);
        data_density->SetName("Density");
        data_velocity->SetNumberOfComponents(3);
        data_velocity->SetName("Velocity");



        for(vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
        {
            double * temp_velocity=&velocity[3*i];
            double * temp_density=&density[i];
            double * temp_phase=&phase[i];
            data_velocity->InsertNextTupleValue(temp_velocity);
            data_density->InsertNextTupleValue(temp_density);
            data_phase->InsertNextTupleValue(temp_phase);
        }

        vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
        structuredGrid->SetDimensions(NX,NY,NZ);
        structuredGrid->SetPoints(points);
        structuredGrid->GetPointData()->AddArray(data_velocity);
        structuredGrid->GetPointData()->AddArray(data_phase);
        structuredGrid->GetPointData()->AddArray(data_density);

        vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
        writer->SetFileName(name.c_str());
        writer->SetInput(structuredGrid);
        writer->Write();
    }

}

void Solver::writeTextXZVelocity(std::string name)
{
    Solver::getWholeVelocity();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterX=0;counterX<NX;counterX++)
                fout << counterX << " " << counterZ << " "<<velocity[3*(counterZ*NX*NY+NY/2*NX+counterX)]<<" "<<
                velocity[3*(counterZ*NX*NY+NY/2*NX+counterX)+1]<<" "<<velocity[3*(counterZ*NX*NY+NY/2*NX+counterX)+2]<<"\n";
    }
}

void Solver::writeTextXYVelocity(std::string name)
{
    Solver::getWholeVelocity();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterY=0;counterY<NY;counterY++)
            for(int counterX=0;counterX<NX;counterX++)
                fout << counterX << " " << counterY << " "<<velocity[3*(NZ/2*NX*NY+counterY*NX+counterX)]<<
                " "<<velocity[3*(NZ/2*NX*NY+counterY*NX+counterX)+1]<<" "<<velocity[3*(NZ/2*NX*NY+counterY*NX+counterX)+2]<<"\n";
    }
}

void Solver::writeTextYZVelocity(std::string name)
{
    Solver::getWholeVelocity();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                fout << counterY << " " << counterZ <<" "<<velocity[3*(counterZ*NX*NY+counterY*NX+NX/2)]<<
                " "<<velocity[3*(counterZ*NX*NY+counterY*NX+NX/2)+1]<<" "<<velocity[3*(counterZ*NX*NY+counterY*NX+NX/2)+2]<<"\n";
    }
}





void Solver::writeTextWholePhase(std::string name)
{
    Solver::getWholePhase();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                for(int counterX=0;counterX<NX;counterX++)
                    fout << counterX << " " << counterY << " " << counterZ << " "<<phase[counterZ*NX*NY+counterY*NX+counterX]<<"\n";
    }
}

//Write planes of phase
void Solver::writeTextXYPhase(std::string name)
{
    Solver::getWholePhase();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterY=0;counterY<NY;counterY++)
            for(int counterX=0;counterX<NX;counterX++)
                fout << counterX << " " << counterY << " "<<phase[NZ/2*NX*NY+counterY*NX+counterX]<<"\n";
    }
}



void Solver::writeTextXZPhase(std::string name)
{
    Solver::getWholePhase();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterX=0;counterX<NX;counterX++)
                fout << counterX << " " << counterZ << " "<<phase[counterZ*NX*NY+NY/2*NX+counterX]<<"\n";
    }
}



void Solver::writeTextYZPhase(std::string name)
{
    Solver::getWholePhase();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                fout << counterY << " " << counterZ <<" "<<phase[counterZ*NX*NY+counterY*NX+NX/2]<<"\n";
    }
}


//Write planes of density
void Solver::writeTextXYDensity(std::string name)
{
    Solver::getWholeDensity();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterY=0;counterY<NY;counterY++)
            for(int counterX=0;counterX<NX;counterX++)
                fout << counterX << " " << counterY << " "<<density[NZ/2*NX*NY+counterY*NX+counterX]<<"\n";
    }
}



void Solver::writeTextXZDensity(std::string name)
{
    Solver::getWholeDensity();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterX=0;counterX<NX;counterX++)
                fout << counterX << " " << counterZ << " "<<density[counterZ*NX*NY+NY/2*NX+counterX]<<"\n";
    }
}



void Solver::writeTextYZDensity(std::string name)
{
    Solver::getWholeDensity();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                fout << counterY << " " << counterZ <<" "<<density[counterZ*NX*NY+counterY*NX+NX/2]<<"\n";
    }
}


//Write velocity
void Solver::writeTextWholeVelocity(std::string name)
{
    Solver::getWholeVelocity();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                for(int counterX=0;counterX<NX;counterX++)
                    fout << counterX << " " << counterY << " " << counterZ << " "<<velocity[3*(counterZ*NX*NY+counterY*NX+counterX)]<<" "
                        <<velocity[3*(counterZ*NX*NY+counterY*NX+counterX)+1]<<" "<<velocity[3*(counterZ*NX*NY+counterY*NX+counterX)+2]<<"\n";
    }
}


void Solver::writeVTKWholePhase(std::string name)
{
    Solver::getWholePhase();
    if (rank==0)
    {
   		name=name+".vtk";
        std::ofstream fout(name.c_str());
        fout<<"# vtk DataFile Version 3.0\n";
        fout<<"Binary liquid phase field vtk representation\n";
        fout<<"ASCII\n\n";
        fout<<"DATASET STRUCTURED_GRID\n";
        fout<<"DIMENSIONS "<<NX<<" "<<NY<<" "<<NZ<<"\n";
        //fout<<"ORIGIN "<<floor(NX/2)<<" "<<floor(NY/2)<<" "<<floor(NZ/2)<<"\n";
        //fout<<"SPACING 1 1 1\n";
        fout<<"POINTS "<<NX*NY*NZ<<" double\n";
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                for(int counterX=0;counterX<NX;counterX++)
                    fout<<counterX<<" "<<counterY<<" "<<counterZ<<"\n";
        fout<<"\n";
        fout<<"POINT_DATA "<<NX*NY*NZ<<"\n";
        fout<<"SCALARS phase double\n";
        fout<<"LOOKUP_TABLE phase_table\n";
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                for(int counterX=0;counterX<NX;counterX++)
                    fout<<phase[counterZ*NX*NY+counterY*NX+counterX]<<"\n";

    }
}



void Solver::writeLocalDensity(std::string name)
{
    std::stringstream process;
    process<<"_process"<<rank;
    name=name+process.str();
    lattice->writeDensity(name);
}
void Solver::writeLocalPhase(std::string name)
{
    std::stringstream process;
    process<<"_process"<<rank;
    name=name+process.str();
    lattice->writePhase(name);
}

void Solver::writeLocalTextDensity(std::string name)
{
    std::stringstream process;
    process<<"_process"<<rank;
    name=name+process.str();
    lattice->writeTextDensity(name);
}
void Solver::writeLocalTextPhase(std::string name)
{
    std::stringstream process;
    process<<"_process"<<rank;
    name=name+process.str();
    lattice->writeTextPhase(name);
}



void Solver::updateMacro()
{
    lattice->updateMacro();
}
void Solver::updatePhase()
{
    //Exchange phase fields to calculate it properly
    int zbegin=lattice->getZbegin();
    int zend=lattice->getZend();
    MPI_Status status;

	for (int iCount=0; iCount<NX*NY; iCount++)
	{
		lattice->phase_top[iCount]=lattice->phase[(zend-zbegin)*NX*NY+iCount];
		lattice->phase_bottom[iCount]=lattice->phase[iCount];
	}

	if (size!=1)
	{
		if(rank==0)
		{
			MPI_Sendrecv_replace(lattice->phase_top,NX*NY,MPI_DOUBLE,1,1,1,1,MPI_COMM_WORLD,&status);
			MPI_Sendrecv_replace(lattice->phase_bottom,NX*NY,MPI_DOUBLE,mpi_wrapper().get_size()-1,1,mpi_wrapper().get_size()-1,1,MPI_COMM_WORLD,&status);
		}
        else
		{
			MPI_Sendrecv_replace(lattice->phase_bottom,NX*NY,MPI_DOUBLE,mpi_wrapper().get_rank()-1,1,mpi_wrapper().get_rank()-1,1,MPI_COMM_WORLD,&status);
			MPI_Sendrecv_replace(lattice->phase_top,NX*NY,MPI_DOUBLE,(mpi_wrapper().get_rank()+1)%mpi_wrapper().get_size(),1,(mpi_wrapper().get_rank()+1)%mpi_wrapper().get_size(),1,MPI_COMM_WORLD,&status);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else
	{
		for (int iCount=0; iCount<NX*NY; iCount++)
		{
			lattice->phase_top[iCount]=lattice->phase[iCount];
			lattice->phase_bottom[iCount]=lattice->phase[(zend-zbegin)*NX*NY+iCount];
		}

	}


}

void Solver::updateWall()
{
	const double phase_wall=params_list("phase_wall").value<double>();
	const double rho_wall=params_list("rho_wall").value<double>();
	const int zbegin=lattice->getZbegin();
	const int zend=lattice->getZend();
	for (int iZ=0;iZ<NZ;iZ++)
		for(int iY=0;iY<NY;iY++)
			for(int iX=0;iX<NX;iX++)
				if ((iZ>=zbegin)&&(iZ<=zend)&&(geom->getType(iX,iY,iZ)==SolidNode))
				{
					putPhase(iX,iY,iZ,iX,iY,iZ,phase_wall);
					putDensity(iX,iY,iZ,iX,iY,iZ,rho_wall);
				}
}

void Solver::init()
{
    updateWall();
	updatePhase();
    lattice->init();
}

void Solver::collide_stream()
{
    updateMacro();
    updatePhase();
    lattice->collide_stream();
    if (size!=1)
    {
        propagate();
    }
    exchangeMatrices();
}

void Solver::propagate()
{
    lattice->preparePopulations();

    MPI_Status status;
	if(mpi_wrapper().get_rank()==0)
	{
		MPI_Sendrecv_replace(lattice->layer_top_f,NX*NY*NPOP,MPI_DOUBLE,1,1,1,1,MPI_COMM_WORLD,&status);
		MPI_Sendrecv_replace(lattice->layer_bottom_f,NX*NY*NPOP,MPI_DOUBLE,mpi_wrapper().get_size()-1,1,mpi_wrapper().get_size()-1,1,MPI_COMM_WORLD,&status);

		MPI_Sendrecv_replace(lattice->layer_top_g,NX*NY*NPOP,MPI_DOUBLE,1,1,1,1,MPI_COMM_WORLD,&status);
		MPI_Sendrecv_replace(lattice->layer_bottom_g,NX*NY*NPOP,MPI_DOUBLE,mpi_wrapper().get_size()-1,1,mpi_wrapper().get_size()-1,1,MPI_COMM_WORLD,&status);
	}
	else
	{
		MPI_Sendrecv_replace(lattice->layer_bottom_f,NX*NY*NPOP,MPI_DOUBLE,mpi_wrapper().get_rank()-1,1,mpi_wrapper().get_rank()-1,1,MPI_COMM_WORLD,&status);
		MPI_Sendrecv_replace(lattice->layer_top_f,NX*NY*NPOP,MPI_DOUBLE,(mpi_wrapper().get_rank()+1)%mpi_wrapper().get_size(),1,(mpi_wrapper().get_rank()+1)%mpi_wrapper().get_size(),1,MPI_COMM_WORLD,&status);

		MPI_Sendrecv_replace(lattice->layer_bottom_g,NX*NY*NPOP,MPI_DOUBLE,mpi_wrapper().get_rank()-1,1,mpi_wrapper().get_rank()-1,1,MPI_COMM_WORLD,&status);
		MPI_Sendrecv_replace(lattice->layer_top_g,NX*NY*NPOP,MPI_DOUBLE,(mpi_wrapper().get_rank()+1)%mpi_wrapper().get_size(),1,(mpi_wrapper().get_rank()+1)%mpi_wrapper().get_size(),1,MPI_COMM_WORLD,&status);
	}

    MPI_Barrier(MPI_COMM_WORLD);

    lattice->finishPropagation();
}

void Solver::exchangeMatrices()
{
    lattice->exchangeMatrices();
}
