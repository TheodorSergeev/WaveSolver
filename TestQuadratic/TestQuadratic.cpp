#include "TestQuadratic.h"

double TestQuadratic::analyt_sol(double x, double t)
{

    return x * (length - x) * (1 + t / 2);

}

double TestQuadratic::InitCoord(double x)
{

    return analyt_sol(x, 0.0);

}

double TestQuadratic::InitVeloc(double x)
{

    return 0.5 * analyt_sol(x, 0.0);

}

double TestQuadratic::HeterFunc(double x, double t)
{

    return 2 * (1 + 0.5 * t) * pow(c, 2);

}

double TestQuadratic::WaveSpeed(double x)
{

    return c;

}

TestQuadratic::TestQuadratic(bool to_print_)
{

    to_print = to_print_;

    cout << "TestQuadratic... ";

    left_bound_cond = DIRICHLET;
    right_bound_cond = DIRICHLET;

    length      = 2.5;
    c           = 1.5;
    courant_num = 0.75;
    mesh_size   = 6;
    time_lim    = 18;
    t_step      = courant_num * length / mesh_size / c;

    x_step      = c * t_step / courant_num;
    mesh_size   = (int) round(length / x_step) + 1;

    curr_sol_arr.assign(mesh_size, 0.0);
    next_sol_arr = curr_sol_arr;
    prev_sol_arr = curr_sol_arr;

    t_st_sq = t_step * t_step;
    coef_sq = t_st_sq / (x_step * x_step);

    Check();
    //Dump();

    Solver();

    cout << "successfull.\n";

}


double TestQuadratic::DirL(double t)
{

    return 0.0;

}

double TestQuadratic::DirR(double t)
{

    return 0.0;

}

