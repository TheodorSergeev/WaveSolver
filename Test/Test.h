#ifndef TEST_H
#define TEST_H

#include "Wave1DSolver.h"
#include "catch2.hpp"

class Test: public Wave1DSolver
{

protected:
    bool to_print;

public:

    virtual double analyt_sol(double x, double t) = 0;
    void PrintAnalyt(double t);
    virtual void user_action(double t);

};

#endif
