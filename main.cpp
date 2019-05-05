#include "Wave1DSolver.h"
#include "Guitar.h"
#include "Test.h"
#include "TestQuadratic.h"
#include "TestConstant.h"

#include <stdio.h>
#include <stdlib.h>
#include <exception>
#include <limits>


// indent

int main()
{

    try
    {

        TestQuadratic test1;
        TestConstant  test2;

        Guitar prob;
        prob.Solver();

    }
    catch(string& err)
    {

        cout << "Exception: " << err << "\n";

    }
    catch(const char* err)
    {

        cout << "Exception: " << err << "\n";

    }
    catch(...)
    {


        cout << "Unknown exception\n" << "\n";

    }

    return 0;

}
