#ifndef TESTQUADRATIC_H
#define TESTQUADRATIC_H

#include "Test.h"

// Check that u(x,t)=x(L-x)(1+t/2) is exactly reproduced.
class TestQuadratic: public Test
{

private:

    double c;

public:

    virtual double analyt_sol(double x, double t);

    virtual double InitCoord(double x);
    virtual double InitVeloc(double x);
    virtual double HeterFunc(double x, double t);
    virtual double WaveSpeed(double x);

    TestQuadratic(bool to_print_ = false);

    virtual double DirL(double t);
    virtual double DirR(double t);

};

#endif
