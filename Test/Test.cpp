#include "Test.h"

void Test::PrintAnalyt(double t)
{

    cout.precision(6);

    cout << "t = " << t << "\n";

    cout << "Analyt:  ";

    for(unsigned int i = 0; i < mesh_size; ++i)
    {

        cout << analyt_sol(i * x_step, t) << " ";

    }

    cout << "\n";
    cout << "Numerical:  ";

    for(unsigned int i = 0; i < mesh_size; ++i)
    {

        cout << curr_sol[i] << " ";

    }

    cout << "\n";

}

void Test::user_action(double t)
{

    if(to_print)
    {

        PrintAnalyt(t);

    }

    for(unsigned int i = 0; i < mesh_size; ++i)
    {

        double anl_sol = analyt_sol(i * x_step, t);
        double num_sol = curr_sol[i];
        const double eps = pow(10, -13);

        //cout << "error: " << fabs(anl_sol - num_sol) << " " << eps << "\n";
        // std::numeric_limits<double>::epsilon() << "\n";

        if(fabs(anl_sol - num_sol) > eps)
        {

            throw "x = " + to_string(i * x_step) +
                  " t = " + to_string(t)          + " :" +
                  " anl_sol = " + to_string(anl_sol) +
                  " !=  num_sol = " + to_string(num_sol);

        }

    }

}
