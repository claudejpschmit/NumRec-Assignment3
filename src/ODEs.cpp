#include "ODEs.hpp"
#include <cmath>

ODE::~ODE()
{
}

double ODE::first_derivative(double y, double x)
{
    (void)y;
    (void)x;
    return 0;
}

double ODE::initial_value()
{
    return 0;
}

double ODE::exact_solution(double x)
{
    (void)x;
    return 0;
}


/* ****      Electric Field     **** */


ElectricFieldODE::ElectricFieldODE(double x0, double h)
    :
    h(h)
{
    this->x0 = x0;
}
double ElectricFieldODE::first_derivative(double y, double x)
{
    (void)y;
    if (x >= 1.0 && x <= 2.0)
        return h;
    else if (x > 2.0 && x <= 3.0)
        return -h;
    else
        return 0.0;
}
    
double ElectricFieldODE::initial_value()
{
    return 0.0;
}

double ElectricFieldODE::exact_solution(double x)
{
    if (x >= 1.0 && x <= 2.0)
        return h * x - h;
    else if (x > 2.0 && x <= 3.0)
        return -h * x + 3 * h;
    else
        return 0.0;
}


/* ****     Potential ODE       **** */

PotentialODE::PotentialODE(double x0, double h, std::vector<double> E_field, double E_field_dx)
    :
    h(h),
    dx(E_field_dx)
{
    this->E_field = E_field;
    this->x0 = x0;
}

double PotentialODE::first_derivative(double y, double x)
{
    (void)y;
    return - E_field[x/dx];
}
double PotentialODE::initial_value()
{
    return 0.0;
}
double PotentialODE::exact_solution(double x)
{
    if (x < 1.0)
        return 0.0;
    else if (x >= 1.0 && x <= 2.0)
        return - 0.5 * h * x * x + h * x - 0.5 * h;
    else if (x > 2.0 && x <= 3.0)
        return 0.5 * h * x * x - 3 * h * x + 3.5 * h;
    else
        return  exact_solution(3);
}
