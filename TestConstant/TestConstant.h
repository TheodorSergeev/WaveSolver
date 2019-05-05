#ifndef TESTCONSTANT_H
#define TESTCONSTANT_H

#include "Test.h"

class TestConstant: public Test
{

private:

    double c;
    double constant_sol;

public:

    virtual double analyt_sol(double x, double t);
    virtual double InitCoord(double x);
    virtual double WaveSpeed(double x);

    TestConstant(bool to_print_ = false);

    virtual double DirL(double t);
    virtual double DirR(double t);

};

#endif
