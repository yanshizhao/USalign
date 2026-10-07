#include "geometry.h"

void transform(const Vec3& t, const RotMat& u, const Vec3& x, Vec3& x1)
{
    x1[0]=t[0]+dot(u[0], x);
    x1[1]=t[1]+dot(u[1], x);
    x1[2]=t[2]+dot(u[2], x);
}

void do_rotation(CoordArray& x, CoordArray& x1, int len, const Vec3& t, const RotMat& u)
{
    for(int i=0; i<len; i++)
    {
        transform(t, u, x[i], x1[i]);
    }
}