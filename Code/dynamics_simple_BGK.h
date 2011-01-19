#ifndef DYNAMICS_SIMPLE_BGK_H
#define DYNAMICS_SIMPLE_BGK_H

#include "dynamics.h"
#include "params_list.h"

template<typename D> class DynamicsSimpleBGK:public Dynamics<D>
{
    public:

    explicit DynamicsSimpleBGK(Solver<D> * _solver,ParamsList _params);
    virtual ~DynamicsSimpleBGK()
    {}

    void collide_stream(int iX,int iY,int iZ);
    void init(int iX,int iY, int iZ);
    virtual void updateFields(int iX,int iY,int iZ);

    public:

    double omega;
    double force_x;
    double force_y;
    double force_z;

};
#endif
