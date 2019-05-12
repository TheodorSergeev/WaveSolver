#include "Wave1DSolver.h"

double PmlVars::abs_coef(int layer)
{

    if(layer < 0)
        throw "Layer l = " + to_string(layer) + " but cannot be <0 in PmlVars::abs_coef";

    return max_abs_coef * pow(layer / layers_num, power);

}

void Wave1DSolver::CheckX(double x)
{

    if(x < 0 || x > length)
        throw "Coordinate x = " + to_string(x) + " but cannot be <0 or > length";

}

double Wave1DSolver::avr_c_sq(double node)
{

    return (pow(WaveSpeed((node + 0.5) * x_step), 2) + pow(WaveSpeed((node - 0.5) * x_step), 2)) / 2;

}

void Wave1DSolver::Check()
{

    if(next_sol.size() != mesh_size)
        throw "mesh_size = " + to_string(mesh_size) + " is not equal to next_sol.size() = " + to_string(next_sol.size());
    if(curr_sol.size() != mesh_size)
        throw "mesh_size = " + to_string(mesh_size) + " is not equal to curr_sol.size() = " + to_string(curr_sol.size());
    if(prev_sol.size() != mesh_size)
        throw "mesh_size = " + to_string(mesh_size) + " is not equal to prev_sol.size() = " + to_string(prev_sol.size());

    if(length      <= 0)
        throw "length = "      + to_string(length)      + " cannot be < 0";
    if(time_lim    <= 0)
        throw "time_lim = "    + to_string(time_lim)    + " cannot be < 0";
    if(t_step      <= 0)
        throw "t_step = "      + to_string(t_step)      + " cannot be < 0";
    if(x_step      <= 0)
        throw "x_step = "      + to_string(x_step)      + " cannot be < 0";
    if(courant_num <= 0)
        throw "courant_num = " + to_string(courant_num) + " cannot be < 0";

    if(left_bound_cond == NONE)
        throw "left_bound_cond must be defined";
    if(right_bound_cond == NONE)
        throw "right_bound_cond must be defined";

}

void Wave1DSolver::Dump()
{

    cout << "length = "      << length      << "\n"
         << "time_lim = "    << time_lim    << "\n"
         << "t_step = "      << t_step      << "\n"
         << "x_step = "      << x_step      << "\n"
         << "courant_num = " << courant_num << "\n"
         << "mesh_size = "   << mesh_size   << "\n";

    cout << "prev_sol:\n";
    for(auto el : prev_sol)
        cout << el << " ";
    cout << "\ncurr_sol:\n";
    for(auto el : curr_sol)
        cout << el << " ";
    cout << "\nnext_sol:\n";
    for(auto el : next_sol)
        cout << el << " ";

    cout << "\n";

}

double Wave1DSolver::LeftBoundCond(double t)
{

    switch(left_bound_cond)
    {

        case NONE:      throw "Left boundary condition must is not defined";
        case DIRICHLET: return DirL(t);
        case NEUMANN:   return NeuL(t);
        case MUR:       return MurL(t);
        case PML:       return PmlL(t);

    }

}

double Wave1DSolver::RightBoundCond(double t)
{

    switch(right_bound_cond)
    {

        case NONE:      throw "Right boundary condition must is not defined";
        case DIRICHLET: return DirR(t);
        case NEUMANN:   return NeuR(t);
        case MUR:       return MurR(t);
        case PML:       return PmlR(t);

    }

}


double Wave1DSolver::DirL(double t)
{

    return 0.0;

}

double Wave1DSolver::DirR(double t)
{

    return 0.0;

}

double Wave1DSolver::NeuL(double t)
{

    double left_bound = 0.0;
    int i = 0;

    if(t == 0.0)
    {

        left_bound = prev_sol[i] + t_step * InitVeloc(0.0) +
                     coef_sq * (avr_c_sq(i + 0.5) * (curr_sol[i + 1] - curr_sol[i]) / 2 -
                                avr_c_sq(i + 0.5) * (curr_sol[i] - curr_sol[i + 1]) / 2) +
                     t_st_sq * HeterFunc(0.0, t) / 2;

    }

    left_bound = - prev_sol[i] + 2 * curr_sol[i] +
                 coef_sq * (avr_c_sq(i + 0.5) * (curr_sol[i + 1] - curr_sol[i]) -
                            avr_c_sq(i + 0.5) * (curr_sol[i] - curr_sol[i + 1])) / 2 +
                 t_st_sq * HeterFunc(0.0, t);

    return left_bound;

}

double Wave1DSolver::NeuR(double t)
{

    double right_bound = 0.0;
    int i = mesh_size - 1;

    if(t == 0.0)
    {

        right_bound = prev_sol[i] + t_step * InitVeloc(x_step * i) +
                      coef_sq * (avr_c_sq(i - 0.5) * (curr_sol[i - 1] - curr_sol[i]) / 2 -
                                 avr_c_sq(i - 0.5) * (curr_sol[i] - curr_sol[i - 1]) / 2) +
                      t_st_sq * HeterFunc(x_step * i, t) / 2;

    }

    right_bound = - prev_sol[i] + 2 * curr_sol[i] +
                  coef_sq * (avr_c_sq(i - 0.5) * (curr_sol[i - 1] - curr_sol[i]) -
                             avr_c_sq(i - 0.5) * (curr_sol[i] - curr_sol[i - 1])) / 2 +
                  t_st_sq * HeterFunc(x_step * i, t);

    return right_bound;

}


double Wave1DSolver::PmlL(double t)
{

    return 0.0;

}

double Wave1DSolver::PmlR(double t)
{

    // todo: check that N < mesh_size / 2 - 1

    for(int i = mesh_size - 1 - pml_right.layers_num; i < mesh_size - 1; ++i)
    {

        // here c = const is considered
        // todo: c = c(x)

        //cout << "i+1 = " << i+1 << " " << mesh_size - 1 << '\n';

        next_sol[i + 1] = curr_sol[i] - t_step / pow(WaveSpeed(i * x_step), 2.0) *
                          (pml_right.abs_coef(mesh_size - 1 - i) * curr_sol[i] + (curr_sol[i] - curr_sol[i - 1]) / x_step);

    }

    return 0.0;

}


double Wave1DSolver::MurL(double t)
{

    //return curr_sol[0] - WaveSpeed(1 * x_step) * t_step / x_step * (curr_sol[1] - curr_sol[0]);

    // the neighbour is already computed
    double ct = WaveSpeed(0.5 * x_step) * t_step;
    return curr_sol[1] - (ct - x_step) / (ct + x_step) * (curr_sol[0] - next_sol[1]);

}

double Wave1DSolver::MurR(double t)
{

    //int N = mesh_size - 1;
    //return curr_sol[N] - WaveSpeed((N - 1) * x_step) * t_step / x_step * (curr_sol[N] - curr_sol[N - 1]);

    // old version - will be revisited
    double ct = WaveSpeed(0.5 * x_step) * t_step;
    return curr_sol[mesh_size - 2] + (ct - x_step) / (ct + x_step) * (next_sol[mesh_size - 2] - curr_sol[mesh_size - 1]);

}

double Wave1DSolver::InitCoord(double x)
{

    CheckX(x);
    return 0.0;

}

double Wave1DSolver::InitVeloc(double x)
{

    CheckX(x);
    return 0.0;

}

double Wave1DSolver::WaveSpeed(double x)
{

    CheckX(x);
    return 0.0;

}

double Wave1DSolver::HeterFunc(double x, double t)
{

    CheckX(x);
    return 0.0;

}

Wave1DSolver::Wave1DSolver():
    length(0), time_lim(0), t_step(0), t_st_sq(0), x_step(0), courant_num(0), coef_sq(0), mesh_size(0),
    next_sol_arr(), curr_sol_arr(), prev_sol_arr(),
    next_sol(next_sol_arr), curr_sol(curr_sol_arr), prev_sol(prev_sol_arr),
    left_bound_cond(NONE), right_bound_cond(NONE)
{}


void Wave1DSolver::RecordCurrSol(string& prefix, double t)
{

    std::ofstream output;
    //string fname = "data/" + prefix + "_" + to_string(t) + ".txt";

    string fname = "data/" + to_string(t) + ".txt";
    output.open(fname, std::ifstream::out);

    if(!output.is_open() || !output.good())
        throw "Bad output filename " + fname + " in RecordCurrSol (check if the directory exists).\n";

    for(unsigned int i = 0; i < mesh_size; ++i)
    {

        output << i * x_step << " " << curr_sol[i] << "\n";

    }

    output.close();

}

void Wave1DSolver::user_action(double t)
{

    RecordCurrSol(pref, t);

}

void Wave1DSolver::Solver()
{

    Check();
    if(mesh_size <= 0) throw "mesh size <=0 in Solver";

    int t_steps_num = (int) round(time_lim / t_step);

    for(unsigned int i = 0; i < mesh_size; ++i)
    {

        prev_sol[i] = InitCoord(x_step * i);
    }

    curr_sol = prev_sol;
    user_action(0.0);

    curr_sol = curr_sol_arr;

    for(int j = 1; j <= t_steps_num; ++j)
    {

        Iteration(j);

    }

}

void Wave1DSolver::Iteration(int t_step_num)
{

    int j = t_step_num;

    if(t_step_num == 1)
    {

        for(unsigned int i = 1; i < mesh_size - 1; ++i)
        {

            curr_sol[i] = prev_sol[i] + t_step * InitVeloc(x_step * i) +
                          coef_sq * (avr_c_sq(i + 0.5) * (prev_sol[i + 1] - prev_sol[i]) -
                                     avr_c_sq(i - 0.5) * (prev_sol[i] - prev_sol[i - 1])) / 2 +
                          t_st_sq * HeterFunc(x_step * i, t_step * (j - 1)) / 2;

        }

        curr_sol[0]             = LeftBoundCond (0.0);
        curr_sol[mesh_size - 1] = RightBoundCond(0.0);

        user_action(j * t_step);

    }
    else
    {

        for(unsigned int i = 1; i < mesh_size - 1; ++i)
        {

            next_sol[i] = - prev_sol[i] +
                          2 * curr_sol[i] +
                          coef_sq * (avr_c_sq(i + 0.5) * (curr_sol[i + 1] - curr_sol[i]) -
                                     avr_c_sq(i - 0.5) * (curr_sol[i] - curr_sol[i - 1])) / 2 +
                          t_st_sq * HeterFunc(x_step * i, t_step * (j - 1)) / 2;

        }

        next_sol[0]             = LeftBoundCond (t_step * (j - 1));
        next_sol[mesh_size - 1] = RightBoundCond(t_step * (j - 1));

        prev_sol = curr_sol;
        curr_sol = next_sol;
        next_sol = prev_sol_arr;

        user_action(j * t_step);

    }

}
