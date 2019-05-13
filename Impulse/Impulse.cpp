#include "Impulse.h"


double Impulse::InitCoord(double x)
{

    //return 0.0;
    CheckX(x);

    if(x <= X0 - (x_step * 10) || x >= X0 + (x_step * 10))
        return 0.0;
    else
        return amp * cos((X0 - x) * M_PI / 2 / (x_step * 10));

}

double Impulse::InitVeloc(double x)
{

    return 0.0;
    CheckX(x);

    if(x < X0 - (x_step * 10) || x > X0 + (x_step * 10))
        return 0.0;
    else
        return - amp * cos((X0 - x) * M_PI * (length * (x_step * 10)));


}

double Impulse::WaveSpeed(double x)
{

    CheckX(x);
    return 1.0;

}

Impulse::Impulse()
{

    length = 0.75;

    double freq = 440;

    double omega = 2 * M_PI * freq;
    periods = 2.0;
    double wavespeed = 1.0;
    amp = 0.005;
    X0 = 0.5 * length;

    t_step      = length / 100.0 / wavespeed;
    courant_num = 0.85;
    x_step      = wavespeed * t_step / courant_num; // !!!
    mesh_size   = (int)(round(length / x_step));

    time_lim    = t_step * 400;

    curr_sol_arr.assign(mesh_size, 0.0);
    next_sol_arr = curr_sol_arr;
    prev_sol_arr = curr_sol_arr;

    t_st_sq = t_step * t_step;
    coef_sq = t_st_sq / (x_step * x_step);

    left_bound_cond  = NEUMANN;
    right_bound_cond = NEUMANN;

    Check();
    //Dump();

}

