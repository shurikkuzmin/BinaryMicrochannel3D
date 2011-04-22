#ifndef DYNAMICS_SIMPLE_SPECIAL
#define DYNAMICS_SIMPLE_SPECIAL
#include <vector>
#include "dynamics_simple_BGK.h"
#include <iostream>
class DynamicsSimpleSpecial : public DynamicsSimpleBGK
{
    public:

    explicit DynamicsSimpleSpecial(Solver * _solver,ParamsList _params);

    virtual ~DynamicsSimpleSpecial()
    {
    }

    void collide_stream(int iX,int iY,int iZ);

    void init(int iX,int iY, int iZ);


    public:

    int iX;
    int iY;
    int iZ;

    std::vector<int> directions;
    int compliment[19];
};

#endif
