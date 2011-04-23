#ifndef DYNAMICS_ZOUHE_SIMPLE_PRESSURE
#define DYNAMICS_ZOUHE_SIMPLE_PRESSURE
#include <vector>
#include "dynamics_BGK.h"
#include <iostream>
class DynamicsZouHeSimplePressure : public DynamicsBGK
{
    public:

    explicit DynamicsZouHeSimplePressure(Solver * _solver,ParamsList _params);

    virtual ~DynamicsZouHeSimplePressure()
    {
    }

    void collide_stream(int iX,int iY,int iZ);

    void init(int iX,int iY, int iZ);
    virtual void updateFields(int iX,int iY,int iZ);

    public:

    int normx;
    int normy;
    int normz;

    double rho_press;
    double phase_press;

    std::vector<int> known;
    std::vector<int> unknown;
    std::vector<int> positive;
    std::vector<int> negative;
    std::vector<int> neutral;
    std::vector<int> diagonal;

    int compliment[19];
};

#endif

