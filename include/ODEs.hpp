#pragma once

#include <iostream>
#include <vector>

class ODE {
        
public:
    virtual ~ODE();
    virtual double first_derivative(double y, double x);
    virtual double initial_value();
    virtual double exact_solution(double x);

protected:
    double x0;
};

class ElectricFieldODE : public ODE {

public:
    ElectricFieldODE(double x0, double h);
    double first_derivative(double y, double x);
    double initial_value();
    double exact_solution(double x);
private:
    double h;
};

class PotentialODE : public ODE {

public:
    PotentialODE(double x0, double h, std::vector<double> E_field, double E_field_dx);
    double first_derivative(double y, double x);
    double initial_value();
    double exact_solution(double x);
    
private:
    double h;
    double dx;
    std::vector<double> E_field;
};

