#include "lattice.h"
#include "dynamics_BGK.h"
#include "mpi_singleton.h"


template<typename D>
void Lattice<D>::writePhase(std::string name)
{

    name=name+".vts";

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    for(int counterZ=0;counterZ<zend-zbegin+1;counterZ++)
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
    structuredGrid->SetDimensions(NX,NY,zend-zbegin+1);
    structuredGrid->SetPoints(points);
    structuredGrid->GetPointData()->AddArray(data);

    vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer->SetFileName(name.c_str());
    writer->SetInput(structuredGrid);
    writer->Write();
}

template<typename D>
void Lattice<D>::writeTextPhase(std::string name)
{

    name=name+".dat";
    std::ofstream fout(name.c_str());
    for(int counterZ=0;counterZ<zend-zbegin+1;counterZ++)
        for(int counterY=0;counterY<NY;counterY++)
            for(int counterX=0;counterX<NX;counterX++)
                fout << counterX << " " << counterY << " " << counterZ << " "<<phase[counterZ*NX*NY+counterY*NX+counterX]<<"\n";

}

template<typename D>
void Lattice<D>::writeTextDensity(std::string name)
{

    name=name+".dat";
    std::ofstream fout(name.c_str());
    for(int counterZ=0;counterZ<zend-zbegin+1;counterZ++)
        for(int counterY=0;counterY<NY;counterY++)
            for(int counterX=0;counterX<NX;counterX++)
                fout << counterX << " " << counterY << " " << counterZ << " "<<rho[counterZ*NX*NY+counterY*NX+counterX]<<"\n";

}

template<typename D>
void Lattice<D>::writeDensity(std::string name)
{

    name=name+".vts";

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    for(int counterZ=0;counterZ<zend-zbegin+1;counterZ++)
        for(int counterY=0;counterY<NY;counterY++)
            for(int counterX=0;counterX<NX;counterX++)
                points->InsertNextPoint(counterX,counterY,counterZ);


    vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
    data->SetNumberOfComponents(1);
    data->SetName("Density");



    for(vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
    {
        double * temp_data=&rho[i];
        data->InsertNextTupleValue(temp_data);
    }

    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
    structuredGrid->SetDimensions(NX,NY,zend-zbegin+1);
    structuredGrid->SetPoints(points);
    structuredGrid->GetPointData()->AddArray(data);

    vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer->SetFileName(name.c_str());
    //writer.SetPoints();
    writer->SetInput(structuredGrid);
    writer->Write();
}

template<typename D>
void Lattice<D>::init()
{
	for (int counter=0;counter<NUMLOCAL;counter++)
	{
		int iZ=counter/(NX*NY);
		int iY=(counter%(NX*NY))/NX;
		int iX=(counter%(NX*NY))%NX;
		dynamics_list[counter]->init(iX,iY,iZ);
	}
}

template<typename D>
void Lattice<D>::collide_stream()
{
	for (int counter=0;counter<NUMLOCAL;counter++)
	{
		int iZ=counter/(NX*NY);
		int iY=(counter%(NX*NY))/NX;
		int iX=(counter%(NX*NY))%NX;
		dynamics_list[counter]->collide_stream(iX,iY,iZ);
	}
}

template<typename D>
void Lattice<D>::preparePopulations()
{
   	for (int iCount=0; iCount<NX*NY; iCount++)
	{
		for (int iPop=0; iPop<D::NPOP; iPop++)
		{
			layer_top_f[iCount*D::NPOP+iPop]=f[NX*NY*(zend-zbegin)*D::NPOP+iCount*D::NPOP+iPop];
			layer_bottom_f[iCount*D::NPOP+iPop]=f[iCount*D::NPOP+iPop];
			layer_top_g[iCount*D::NPOP+iPop]=g[NX*NY*(zend-zbegin)*D::NPOP+iCount*D::NPOP+iPop];
			layer_bottom_g[iCount*D::NPOP+iPop]=g[iCount*D::NPOP+iPop];
		}
	}
}

template<typename D>
void Lattice<D>::finishPropagation()
{
   	for (int iCount=0; iCount<NX*NY; iCount++)
	{
        int iY=iCount/NX;
        int iX=iCount%NX;
        for(int k=0;k<D::NPOP;k++)
        {
			int iY2=(iY+D::cy[k]+NY)%NY;
			int iX2=(iX+D::cx[k]+NX)%NX;
            if (geometry[(zend-zbegin)*NX*NY+iY2*NX+iX2]==FluidNode)
            {
                //cout<<mpi_wrapper().get_rank()<<"\n";
                if ((D::cz[k]<0) && (geom_top[iY*NX+iX]==FluidNode))
                {
                    f2[((zend-zbegin)*NX*NY+iY2*NX+iX2)*D::NPOP+k]=layer_top_f[(iY*NX+iX)*D::NPOP+k];
                    g2[((zend-zbegin)*NX*NY+iY2*NX+iX2)*D::NPOP+k]=layer_top_g[(iY*NX+iX)*D::NPOP+k];
                }
            }
            if (geometry[iY2*NX+iX2]==FluidNode)
            {
                if ((D::cz[k]>0) && (geom_bottom[iY*NX+iX]==FluidNode) )
                {
                    f2[(iY2*NX+iX2)*D::NPOP+k]=layer_bottom_f[(iY*NX+iX)*D::NPOP+k];
                    g2[(iY2*NX+iX2)*D::NPOP+k]=layer_bottom_g[(iY*NX+iX)*D::NPOP+k];
                }
            }
        }

	}
}

template<typename D>
void Lattice<D>::exchangeMatrices()
{
		//cout<<"I am here\n";
		//char ch;
		//cin>>ch;
		double* pointer;
        //cout<<"I exchange those matrixes"<<f<<" and "<<f2<<"\n";

		pointer=f;
		f=f2;
		f2=pointer;

		pointer=g;
		g=g2;
		g2=pointer;

}
