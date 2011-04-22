#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "params_list.h"
//#include "solver.h"
class Solver;

class Dynamics
{
    public:

    explicit Dynamics(Solver * _solver,ParamsList _params):
        solver(_solver),params(_params)
    {}

    virtual ~Dynamics()
    {}

    virtual void collide_stream(int iX,int iY,int iZ)=0;
	virtual void init(int iX,int iY,int iZ)=0;

    virtual void updateFields(int iX,int iY,int iZ)=0;

    public:
    Solver * solver;
    ParamsList params;

};
#endif
