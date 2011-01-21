#include "dynamics_symmetric.h"
#include "solver.h"
template<typename D>
DynamicsSymmetric<D>::DynamicsSymmetric(Solver<D> * _solver, ParamsList _params_list):
    DynamicsBGK<D>(_solver,_params_list)
{

   	int NX=this->solver->geom->getNX();
	int NY=this->solver->geom->getNY();
    int NZ=this->solver->geom->getNZ();

    coorX=this->solver->iterX;
    coorY=this->solver->iterY;
    coorZ=this->solver->iterZ;

    sourceZ=coorZ;

    bool x=false;
    bool y=false;

    if (coorX+1==NX)
    {
        x=true;
        sourceX=NX-2;
        sourceY=coorY;
    }
    else if(coorX-1==-1)
    {
        x=true;
        sourceX=1;
        sourceY=coorY;
    }
    else if(coorY+1==NY)
    {
        y=true;
        sourceX=coorX;
        sourceY=NY-2;
    }
    else if(coorY-1==-1)
    {
        y=true;
        sourceX=coorX;
        sourceY=1;
    }

    symmetric[0]=0;
    if(x==true)
        for(int k=1;k<D::NPOP-1;k++)
            for(int l=1;l<D::NPOP-1;l++)
                if ((-D::cx[k]==D::cx[l])&&(D::cy[k]==D::cy[l])&&(D::cz[k]==D::cz[l]))
                {
                    symmetric[k]=l;
                    break;
                }

    if(y==true)
        for(int k=1;k<D::NPOP-1;k++)
            for(int l=1;l<D::NPOP-1;l++)
                if ((D::cx[k]==D::cx[l])&&(D::cy[k]==-D::cy[l])&&(D::cz[k]==D::cz[l]))
                {
                    symmetric[k]=l;
                    break;
                }

    if (x&&y)
    {
        for(int k=1;k<D::NPOP-1;k++)
            for(int l=1;l<D::NPOP-1;l++)
                if ((D::cx[k]==-D::cx[l])&&(D::cy[k]==-D::cy[l])&&(D::cz[k]==D::cz[l]))
                {
                    symmetric[k]=l;
                    break;
                }
        if (coorX==NX-2)
        {
            sourceX=NX-2;
            sourceY=NY-2;
        }
        else
        {
            sourceX=1;
            sourceY=1;
        }
    }

    for(int k=0;k<D::NPOP;k++)
        std::cout<<"Symmetric[k]="<<symmetric[k]<<"\n";

    bounce_back=false;

    Geometry * geom=this->solver->geom;
    for(int k=0; k<D::NPOP; k++)
    {
        int iX2=(sourceX+D::cx[k]+NX)%NX;
        int iY2=(sourceY+D::cy[k]+NY)%NY;
        int iZ2=(sourceZ+D::cz[k]+NZ)%NZ;
        if (geom->getType(iX2,iY2,iZ2)==SolidNode)
            bounce_back=true;
    }

	Lattice<D> * lattice=this->solver->getLattice();

    if (bounce_back)
        directions=lattice->dynamics_list[sourceZ*NX*NY+sourceY*NX+sourceX]->getDirections();

}

template<typename D>
void DynamicsSymmetric<D>::init(int iX,int iY,int iZ)
{

	Lattice<D> * lattice=this->solver->getLattice();
    int zbegin=lattice->getZbegin();
    int zend=lattice->getZend();
    int NX=this->solver->geom->getNX();
    int NY=this->solver->geom->getNY();

    int counter=iZ*NX*NY+iY*NX+iX;
    int source=sourceZ*NX*NY+sourceY*NX+sourceX;


    for(int k=0;k<D::NPOP;k++)
    {
        lattice->f[counter*D::NPOP+k]=lattice->f[source*D::NPOP+symmetric[k]];
		lattice->g[counter*D::NPOP+k]=lattice->g[source*D::NPOP+symmetric[k]];
    }

}



template<typename D>
void DynamicsSymmetric<D>::collide_stream(int iX,int iY, int iZ)
{

	Lattice<D> * lattice=this->solver->getLattice();
    int zbegin=lattice->getZbegin();
    int zend=lattice->getZend();
    int NX=this->solver->geom->getNX();
    int NY=this->solver->geom->getNY();

    int counter=iZ*NX*NY+iY*NX+iX;
    int source=sourceZ*NX*NY+sourceY*NX+sourceX;


    for(int k=0;k<D::NPOP;k++)
    {
        lattice->f[counter*D::NPOP+k]=lattice->f[source*D::NPOP+symmetric[k]];
		lattice->g[counter*D::NPOP+k]=lattice->g[source*D::NPOP+symmetric[k]];
    }

	//Streaming procedure
	for (int k=0; k<D::NPOP; k++)
	{
   		bool flag=false;
   		for (int j=0;j<directions.size();j++)
   		{
            if(k==directions[j])
            {
                flag=true;
				//cout<<"BB\n";
                lattice->f2[counter*D::NPOP+D::compliment[k]]=lattice->f[counter*D::NPOP+k];
                lattice->g2[counter*D::NPOP+D::compliment[k]]=lattice->g[counter*D::NPOP+k];

            }
        }
   		if (!flag)
   		{
            int iZ2=(iZ+D::cz[k]+zend-zbegin+1)%(zend-zbegin+1);
            int iY2=(iY+D::cy[k]+NY)%NY;
            int iX2=(iX+D::cx[k]+NX)%NX;
            if (this->solver->geom->getType(iX,iY,iZ)!=SolidNode)
            {
                lattice->f2[(iZ2*NX*NY+iY2*NX+iX2)*D::NPOP+k]=lattice->f[counter*D::NPOP+k];
                lattice->g2[(iZ2*NX*NY+iY2*NX+iX2)*D::NPOP+k]=lattice->g[counter*D::NPOP+k];
            }
   		}
    }

}

