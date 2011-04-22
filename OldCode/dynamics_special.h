#ifndef DYNAMICS_SPECIAL
#define DYNAMICS_SPECIAL
#include <vector>
#include "dynamics_BGK.h"
#include <iostream>
class DynamicsSpecial : public DynamicsBGK
{
    public:

    explicit DynamicsSpecial(Solver * _solver,ParamsList _params);

    virtual ~DynamicsSpecial()
    {
    }

    void collide_stream(int iX,int iY,int iZ);

    void init(int iX,int iY, int iZ);

    void update_wall_phase(int iX,int iY,int iZ);
    void update_all_gradients(int iX,int iY,int iZ);


    public:

    int iX;
    int iY;
    int iZ;

    double normx;
    double normy;
    double normz;

    double grad[7];

    double phase_gradient;

    std::vector<int> directions;
    std::vector<double> wall_densities;
    int compliment[19];
};

#endif
