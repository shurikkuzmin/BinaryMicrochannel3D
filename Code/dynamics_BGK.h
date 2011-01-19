#ifndef DYNAMICS_BGK_H
#define DYNAMICS_BGK_H

#include "dynamics.h"
#include "params_list.h"

template<typename D>
class DynamicsBGK:public Dynamics<D>
{
    public:

    explicit DynamicsBGK(Solver<D> * _solver,ParamsList _params);
    virtual ~DynamicsBGK()
    {}

    void collide_stream(int iX,int iY,int iZ);
    void init(int iX,int iY, int iZ);
    virtual void updateFields(int iX,int iY,int iZ);

    public:

    double aconst;
    double kconst;
    double gammaconst;
    double tau_phi;
    double tau_liq;
    double tau_gas;
    double force_x;
    double force_y;
    double force_z;

};
#endif
