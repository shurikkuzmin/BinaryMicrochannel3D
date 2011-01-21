#ifndef DYNAMICS_SYMMETRIC
#define DYNAMICS_SYMMETRIC
#include <vector>
#include "dynamics_BGK.h"
#include <iostream>
template<typename D>
class DynamicsSymmetric : public DynamicsBGK<D>
{
    public:

    explicit DynamicsSymmetric(Solver<D> * _solver,ParamsList _params);

    virtual ~DynamicsSymmetric()
    {
    }

    void collide_stream(int iX,int iY,int iZ);

    void init(int iX,int iY, int iZ);

    void update_wall_phase(int iX,int iY,int iZ);
    void update_all_gradients(int iX,int iY,int iZ);


    public:

    int sourceX;
    int sourceY;
    int sourceZ;

    int coorX;
    int coorY;
    int coorZ;

    bool bounce_back;

    std::vector<int> directions;
    int symmetric[D::NPOP];
};

#endif
