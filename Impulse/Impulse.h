#ifndef IMPULSE_H
#define IMPULSE_H

#include "Wave1DSolver.h"

class Impulse: public Wave1DSolver
{

private:

    double amp, periods, X0, a;

public:

    virtual double WaveSpeed(double x);
    virtual double InitCoord(double x);

    Impulse();

};

#endif

