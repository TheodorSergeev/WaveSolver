#ifndef WAVE1DSOLVER_H
#define WAVE1DSOLVER_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using std::vector;
using std::string;
using std::cout;
using std::to_string;

enum BoundCondType { NONE, DIRICHLET, NEUMANN, MUR, PML };

struct PmlVars
{

    int    layers_num;   // number of PML Layers
    double max_abs_coef; // constant determining spatial distribution pattern of abs_coef
    double power;        // maximum possible value of abs_coef

    double abs_coef(int layer); // absorbing coefficient

};

class Wave1DSolver
{

protected:
    string pref;

    vector <double> next_sol_arr; // t+1
    vector <double> curr_sol_arr; // t
    vector <double> prev_sol_arr; // t-1

    vector <double>& next_sol;    // references to arrays above
    vector <double>& curr_sol;    // for avoiding copying when swaping arrays
    vector <double>& prev_sol;

    unsigned int mesh_size;

    BoundCondType left_bound_cond, right_bound_cond;

    double length;
    double time_lim;
    double t_step;
    double courant_num;
    double x_step;

    double t_st_sq, coef_sq;      // utility variables
    double avr_c_sq(double node); // utility functions

    PmlVars pml_left, pml_right;

    void CheckX(double x);        // check whether x is in [0, length]
    void Check();                 // check whether variables are defined correctly
    void Dump();                  // print techical info (all variables)

    double LeftBoundCond (double t);
    double RightBoundCond(double t);

    double NeuL(double t);        // compute left  node solution for Neumann boundary condition
    double NeuR(double t);        // compute right node solution for Neumann boundary condition
    double PmlL(double t);        // compute left  node solution for PML boundary condition
    double PmlR(double t);        // compute right node solution for PML boundary condition
    double MurL(double t);        // compute left  node solution for Mur boundary condition
    double MurR(double t);        // compute left  node solution for Mur boundary condition

public:

    virtual double DirL(double t);
    virtual double DirR(double t);

    virtual double InitCoord(double x);
    virtual double InitVeloc(double x);
    virtual double WaveSpeed(double x);
    virtual double HeterFunc(double x, double t);

    Wave1DSolver();

    void RecordCurrSol(string& prefix, double t);

    virtual void user_action(double t);

    void Solver();
    void Iteration(int t_step_num);

};

#endif
