#ifndef GEOMETRY_H
#define GEOMETRY_H

enum NodeType {FluidNode, SolidNode, PressureNode};

class Geometry
{
    public:

    explicit Geometry(int _NX,int _NY,int _NZ):
        NX(_NX),NY(_NY),NZ(_NZ)
    {
        geometry=new NodeType[NX*NY*NZ];
    }
    ~Geometry()
    {
        delete[] geometry;
    }

    void setType(NodeType type)
    {
        setType(0,0,0,NX-1,NY-1,NZ-1,type);
    }

    void setType(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,NodeType type)
    {
        for(int iZ=iZbegin;iZ<=iZend;iZ++)
            for(int iY=iYbegin;iY<=iYend;iY++)
                for(int iX=iXbegin;iX<=iXend;iX++)
                    geometry[iZ*NX*NY+iY*NX+iX]=type;
    }

    int getNX()
    {
        return NX;
    }
    int getNY()
    {
        return NY;
    }
    int getNZ()
    {
        return NZ;
    }

    NodeType getType(int iX,int iY,int iZ)
    {
         return geometry[iZ*NX*NY+iY*NX+iX];
    }
    private:

    NodeType* geometry;
    int NX;
    int NY;
    int NZ;
};

#endif
