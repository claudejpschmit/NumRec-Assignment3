#pragma once

#include "ODEs.hpp"
#include <iostream>

namespace PNJUNCTION {

    /** \brief A scaffold class to build Ode Solver Methods
     *
     *  ODE Solver Methods should be derived from this class.
     *
     */

    class StepEngine {

        public:
            /** Virtual destructor is necessary to be able 
             *  to delete pointers to derived classes when
             *  they go out of scope
             */
            virtual ~StepEngine();

            /** \brief Contains the ODE solving algorithm. 
             *      To be overloaded by derived classes.
             */
            virtual double step();

        protected:
            /** \brief Pointer to the ODE that is to be solved 
             *      by any of the derived methods.
             */
            ODE* ode; 
            /// y(t_{n-1}) stores the result of the last iteration
            double y_n;
            /// \Delta t is the stepsize between iterations
            double dt;
            /// t_n gives the current time
            double tn; 

    };

    /** \brief Euler Method 
     *
     *  Simple first order method to solve ODEs using interpolation. 
     *  The slope at each point is taken to optain the result.
     */

    class EulerMethod : public StepEngine {

        public:
            /** \brief Constructor for an object to solve a given ODE via
             *      the Euler Method.
             *  \param dt step size for the method
             *  \param ode ODE that is to be solved by the Method
             *  \param t0 initial time
             */
            EulerMethod(double dt, ODE* ode, double t0);

            double step();

    };

    /** \brief Midpoint Runge-Kutta Method
     *
     *  Second order method to solve ODEs. 
     *  Its algorithm is based on the slope at the midpoint of a step interval.
     *  
     */

    class MRKMethod : public StepEngine {

        public:
            /** \brief Constructor for an object to solve a given ODE via
             *      the Midpoint Runge-Kutta Method.
             *  \param dt step size for the method
             *  \param ode ODE that is to be solved by the Method
             *  \param t0 initial time
             */
            MRKMethod(double dt, ODE* ode, double t0);

            double step();

    };

    /** \brief Fourth Order Runge-Kutta Method 
     *
     *  Fourth Order method to solve ODEs. 
     *  Its algorithm is based on the slope at the beginning of the step inteval,
     *  two slope estimates at the midpoint and the slope at the end of this interval.
     *
     */

    class FORKMethod : public StepEngine {

        public:
            /** \brief Constructor for an object to solve a given ODE via
             *      the Fourth Order Runge-Kutta Method.
             *  \param dt step size for the method
             *  \param ode ODE that is to be solved by the Method
             *  \param t0 initial time
             */
            FORKMethod(double dt, ODE* ode, double t0);

            double step();

    };

}
