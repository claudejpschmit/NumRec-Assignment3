#include "ODEs.hpp"
#include "ODE_Solver.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <iomanip>
#include <fftw3.h>


typedef unsigned int uint ;

#define INPUTFORM (1.0 + cos((period/N)*ii)  + cos((3*period/N)*ii)  )
#define XMAX 4.0
#define ITERATIONS 1000
#define DX XMAX/(double)ITERATIONS
#define X0 0.0
#define H 5
#define PI 3.14159265

using namespace std;

int main( ) {

    PNJUNCTION::ElectricFieldODE E_ODE(X0, H);

    PNJUNCTION::FORKMethod F_method_efield(DX, &E_ODE, X0);
    ofstream solution("solution_electric.txt");
    solution << X0 << " " << E_ODE.initial_value() << endl;

    std::vector<double> efield;
    efield.push_back(E_ODE.initial_value());
    for (int i = 1; i < ITERATIONS; ++i) {
        double buffer = F_method_efield.step();
        solution << i * DX << " " << buffer << endl;
        efield.push_back(buffer);
    }
    solution.close();

    PNJUNCTION::PotentialODE P_ODE(X0, H, efield, DX);
    PNJUNCTION::FORKMethod F_method_potential(DX, &P_ODE, X0);

    ofstream solutionP("solution_potential.txt");
    solutionP << X0 << " " << P_ODE.initial_value() << endl;

    for (int i = 1; i < ITERATIONS; ++i) {
        solutionP << i * DX << " " << F_method_potential.step() << endl;
    }
    solutionP.close();

    ofstream solutionE("solution_potential_exact.txt");

    for (int i = 0; i < ITERATIONS; ++i) {
        solutionE << i * DX << " " << P_ODE.exact_solution(i * DX) << endl;
    }
    solutionE.close();
    ofstream solutionEE("solution_electric_exact.txt");

    for (int i = 0; i < ITERATIONS; ++i) {
        solutionEE << i * DX << " " << E_ODE.exact_solution(i * DX) << endl;
    }
    solutionEE.close();


    //////////////////////////////////////////////////// FFTW
    //N = Number of Bins -> more bins = better result
    //
    //TODO: COMPLEX VALUES STILL WRONG
    int N=500 ;
    fftw_complex *rho, *rho_ft, *e_ft, *e;
    fftw_plan p;
    //
    double twopi = 2.0*PI ;
    double period = twopi;
    // this was: double period = twopi * 4.0;
    //....	

    ofstream fft("fft_result.txt");
    rho = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    rho_ft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    e = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    e_ft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    p = fftw_plan_dft_1d(N, rho, rho_ft, FFTW_FORWARD, FFTW_ESTIMATE);

    cout << "Input - REAL - IMAG" << endl;
    for(int n=0; n < N ; ++n) { 
        rho[n][0] = E_ODE.first_derivative(0.0, n * XMAX/(double)N);
        rho[n][1] = 0.0;
        cout << n * (XMAX/(double)N) << " " << rho[n][0] << " " << rho[n][1] << endl;
    }
    //...

    fftw_execute(p); /* repeat as needed */

    //...

    cout << " Results of fftw run " << endl ;
    for(int n=0; n < N ; ++n) { 
        cout << fixed ;
        cout << setprecision(1) ;
        cout << rho_ft[n][0]/(N/2.0) << " " << rho_ft[n][1]/(N/2.0) << endl ;
        // divide rho_ft by ik_n
        e_ft[n][0] = rho_ft[n][1]/(PI * n);
        e_ft[n][1] = - rho_ft[n][0] / (PI * n);
        //os << ii+1 << " "  << out[ii][0]/1000. << endl ;
    }
    e_ft[0][0] = 0;
    e_ft[0][1] = 0;

    p = fftw_plan_dft_1d(N, e_ft, e, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(p); /* repeat as needed */

    cout << " Results of inverse " << endl ;
    for(int n=0; n < N ; ++n) { 
        cout << fixed ;
        cout << setprecision(1) ;
        cout << e[n][0]/N <<  "   "  << e[n][1]/N << endl;
        fft << n * (XMAX/(double)N) << " " << (H/(((e[N/2][0]/(double)N) - (e[0][0]/(double)N)))) * ((e[n][0]/(double)N) - (e[0][0]/(double)N)) << endl;
    }
    fft.close(); 

    fftw_destroy_plan(p);
    fftw_free(rho); fftw_free(rho_ft);


    return 0;


}
