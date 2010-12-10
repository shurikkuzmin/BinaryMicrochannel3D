#ifndef DYNAMICS_BGK_H
#define DYNAMICS_BGK_H

#include "dynamics.h"
#include "params_list.h"

class DynamicsBGK:public Dynamics
{
    public:

    explicit DynamicsBGK(Solver * _solver,ParamsList _params);
    virtual ~DynamicsBGK()
    {}

    void collide_stream(int iX,int iY,int iZ);
    void init(int iX,int iY, int iZ);
    virtual void updateFields(int iX,int iY,int iZ);

    public:

    double omega;
    double aconst;
    double kconst;
    double gammaconst;
    double tau_phi;
    double tau_liq;
    double tau_gas;
    double force_x;
    double force_y;
    double force_z;

    //Constants for the calculations
    double wxx[19];
    double wyy[19];
    double wzz[19];
    double wxy[19];
    double wyz[19];
    double wzx[19];
    double weights[19];

    char cx[19];
    char cy[19];
    char cz[19];
    int NPOP;

    //Parameters for stencils
    double laplace_stencil[19];
    double gradx_stencil[19];
    double grady_stencil[19];
    double gradz_stencil[19];

};
#endif
