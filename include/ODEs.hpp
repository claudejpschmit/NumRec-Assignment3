#pragma once

#include <iostream>
#include <vector>

namespace PNJUNCTION {

    /** \brief A scaffold class to build ODEs from. 
     * 
     *  Functions to be used by an ODE_Solver StepEngine
     *  should derive from this and their definition, 
     *  initial value and exact solution should be
     *  specified. If the exact solution is unknown, 
     *  the function exact_solution(double x) should 
     *  be defined as an empty or zero function. 
     */
    class ODE {

        public:
            /** Virtual destructor is necessary to be able 
             *  to delete pointers to derived classes when 
             *  they go out of scope
             */
            virtual ~ODE();

            /** \brief Contains the definition of the ODE.
             *      It defines what the first derivative of the 
             *      quantity to be solved for is equal to.
             *      This is to be overloaded by the derived classes.
             *  \param y Value of the y parameter of the ODE.
             *  \param x Value of the x parameter of the ODE.
             */
            virtual double first_derivative(double y, double x);

            /** \brief Contains the Initial value of the problem.
             *      To be overloaded by derived classes.
             */
            virtual double initial_value();

            /** \brief Contains the exact solution of the problem if known.
             *      To be overloaded by derived classes. 
             *  \param x Value of the independent parameter in the solution.
             */
            virtual double exact_solution(double x);

        protected:

            /// \brief Stores the starting point for the algorithm.
            double x0;
    };

    /// \brief Electric field equation
    class ElectricFieldODE : public ODE {

        public:

            /** \brief Constructor for Electric field equation
             *  \param x0 position for boundary condition
             *  \param h value of the charge density 
             */
            ElectricFieldODE(double x0, double h);

            /** \brief defines dE/dx = h in [1,2], dE/dx = -h in (2,3] and 0 otherwise.
             *  \param y unused parameter, necessary to overload virtual parent method.
             *  \param x x-coordinate of the p-n junction.
             */
            double first_derivative(double y, double x);

            /** \brief defines boundary condition E(x=0) = 0.
            */
            double initial_value();

            /** \brief defines the analytic solution to the problem, 
             *      used for comparison.
             *  \param x free parameter of solution.
             */
            double exact_solution(double x);
        private:
            /// \brief Local value for charge density.
            double h;
    };


    /// \brief Potential equation
    class PotentialODE : public ODE {

        public:
            /** \brief Constructor for Potential equation.
             *  \param x0 position for the boundary condition.
             *  \param h value of the charge density.
             *  \param E_field vector of values of the electric field,
             *      calculated in the first part of the question.
             *  \param E_field_dx stepsize used to compute E_field values.
             */
            PotentialODE(double x0, double h, std::vector<double> E_field, double E_field_dx);

            /** \brief defines dV/dx = - E(x).
             *  \param y unused parameter, necessary to overload virtual parent method.
             *  \param x x-coordinate of the p-n junction.
             */
            double first_derivative(double y, double x);

            /// \brief defines boundary condition V(x=0) = 0.
            double initial_value();

            /** \brief defines the analytic solution to the problem, 
             *      used for comparison.
             *  \param x free parameter of solution.
             */
            double exact_solution(double x);

        private:
            /// \brief local value for charge density
            double h;
            /// \brief stepsize used to compute Electric field
            double dx;
            /// \brief electric field values
            std::vector<double> E_field;
    };

}
