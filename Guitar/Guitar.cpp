#include "Guitar.h"


double Guitar::InitCoord(double x)
{

    //return 0.0;
    CheckX(x);
    return (x < X0) ? a * x / X0 : a / (length - X0) * (length - x);

}

double Guitar::WaveSpeed(double x)
{

    CheckX(x);
    //if(x > length * 0.1 && x < length * 0.4)
    //    return wavespeed * 2.0;

    return wavespeed;

}

Guitar::Guitar()
{

    length = 0.75;
    freq = 440;
    wavelength = 2 * length;
    wavespeed = freq * wavelength;
    omega = 2 * M_PI * freq;
    periods = 4.0;

    time_lim    = 2 * M_PI / omega * periods;
    t_step      = length / 100.0 / wavespeed;
    courant_num = 0.85;
    x_step      = wavespeed * t_step / courant_num; // !!!

    a = 0.005;
    X0 = 0.8 * length;
    mesh_size = (int)(round(length / x_step));

    curr_sol_arr.assign(mesh_size, 0.0);
    next_sol_arr = curr_sol_arr;
    prev_sol_arr = curr_sol_arr;

    t_st_sq = t_step * t_step;
    coef_sq = t_st_sq / (x_step * x_step);

    left_bound_cond  = DIRICHLET;
    right_bound_cond = DIRICHLET;

    Check();
    //Dump();

}

double Guitar::DirL(double t)
{

    return 0.0;

}

double Guitar::DirR(double t)
{

    return 0.0;

}


