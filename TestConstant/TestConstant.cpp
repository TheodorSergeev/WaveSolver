#include "TestConstant.h"

double TestConstant::analyt_sol(double x, double t)
{

    return constant_sol;

}

double TestConstant::InitCoord(double x)
{

    return analyt_sol(x, 0.0);

}

double TestConstant::WaveSpeed(double x)
{

    return c;

}


TestConstant::TestConstant(bool to_print_)
{

    to_print = to_print_;

    cout << "TestConstant... ";

    constant_sol = 123;

    length      = 2.5;
    c           = 1.5;
    courant_num = 0.75;
    mesh_size   = 6;
    time_lim    = 5;
    t_step      = courant_num * length / mesh_size / c;

    x_step      = c * t_step / courant_num;
    mesh_size   = (int) round(length / x_step) + 1;

    curr_sol_arr.assign(mesh_size, 0.0);
    next_sol_arr = curr_sol_arr;
    prev_sol_arr = curr_sol_arr;

    t_st_sq = t_step * t_step;
    coef_sq = t_st_sq / (x_step * x_step);

    //Dump();

    if(to_print)
        cout << "NEUMANN + NEUMANN \n";

    left_bound_cond  = NEUMANN;
    right_bound_cond = NEUMANN;
    Check();
    Solver();

    if(to_print)
        cout << "DIRICHLET + NEUMANN \n";

    left_bound_cond  = DIRICHLET;
    right_bound_cond = NEUMANN;
    Check();
    Solver();

    if(to_print)
        cout << "NEUMANN + DIRICHLET \n";
    left_bound_cond  = NEUMANN;
    right_bound_cond = DIRICHLET;
    Check();
    Solver();

    if(to_print)
        cout << "DIRICHLET + DIRICHLET \n";
    left_bound_cond  = DIRICHLET;
    right_bound_cond = DIRICHLET;
    Check();
    Solver();

    cout << "successfull.\n";

}

double TestConstant::DirL(double t)
{

    return constant_sol;

}

double TestConstant::DirR(double t)
{

    return constant_sol;

}
