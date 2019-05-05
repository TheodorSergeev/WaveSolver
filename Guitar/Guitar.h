#ifndef GUITAR_H
#define GUITAR_H

#include "Wave1DSolver.h"

class Guitar: public Wave1DSolver
{

private:

    double freq, wavelength, wavespeed, omega, periods, X0, a;

public:

    virtual double InitCoord(double x);
    virtual double WaveSpeed(double x);

    Guitar();

    virtual double DirL(double t);
    virtual double DirR(double t);

};

#endif
