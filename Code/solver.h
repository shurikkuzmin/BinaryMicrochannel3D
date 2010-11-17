#ifndef SOLVER_H
#define SOLVER_H

#include "geometry.h"
#include "lattice.h"
#include "mpi_singleton.h"
#include "dynamics_BGK.h"

#include "destroyer.h"
#include "params_list.h"

class Solver
{
    public:
    explicit Solver(Geometry * _geom,ParamsList _params_list);
    ~Solver()
    {
        delete lattice;
        if (rank==0)
        {
            delete[] density;
            delete[] phase;
            delete[] velocity;
        }
    }

    bool checkNAN();

    void putDensity(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,double dense);
    void putDensity(double dense);
    void putPhase(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,double dense);
    void putPhase(double dense);
    void putVelocity(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,double * velocity);
    void putVelocity(double * velocity);

    void getWholeDensity();
    void getWholePhase();
    void getWholeVelocity();

    void writeWholeDensity(std::string name);
    void writeWholePhase(std::string name);
    void writeWholeVelocity(std::string name);
    void writeWholeDensityPhaseVelocity(std::string name);

    void writeTextXZVelocity(std::string name);
    void writeTextXYVelocity(std::string name);
    void writeTextYZVelocity(std::string name);


    void writeLocalDensity(std::string name);
    void writeLocalPhase(std::string name);
    void writeLocalTextDensity(std::string name);
    void writeLocalTextPhase(std::string name);

    void writeTextWholeVelocity(std::string name);
    void writeTextWholePhase(std::string name);
    void writeVTKWholePhase(std::string name);

    void writeTextXZPhase(std::string name);
    void writeTextXYPhase(std::string name);
    void writeTextYZPhase(std::string name);

    void writeTextXZDensity(std::string name);
    void writeTextXYDensity(std::string name);
    void writeTextYZDensity(std::string name);

    void init();
    void updateWall();
	void updatePhase();
    void updateMacro();
    void propagate();
    void collide_stream();
    void exchangeMatrices();


    Lattice * getLattice()
    {
        return lattice;
    }


	Geometry * geom;

    //Incremental coefficients
    int iterX;
    int iterY;
    int iterZ;

    //Velocity set. TODO: put it in the Descriptor Object
    static const char cx[19];
    static const char cy[19];
    static const char cz[19];

	private:

    const int size;
    const int rank;
    const int NX;
    const int NY;
    const int NZ;
    const int NPOP;


    Lattice * lattice;
    double * density;
    double * phase;
    double * velocity;

    Destroyer<Dynamics> destroyer_list;
    ParamsList params_list;
};

#endif
