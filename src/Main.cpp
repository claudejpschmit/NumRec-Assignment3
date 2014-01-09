#include "ODEs.hpp"
#include "ODE_Solver.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <fftw3.h>

//Convenient variable definitions
#define XMAX 4.0
#define ITERATIONS 1000
#define DX XMAX/(double)ITERATIONS
#define X0 0.0
#define H 5
#define PI 3.14159265

using namespace std;

int main( ) {

    //instantiate electric field ODE and Method to solve it
    PNJUNCTION::ElectricFieldODE E_ODE(X0, H);
    PNJUNCTION::FORKMethod F_method_efield(DX, &E_ODE, X0);

    //creating output file
    ofstream solution("solution_electric.txt");
    solution << X0 << " " << E_ODE.initial_value() << endl;

    //creating vector to contain values of the electric field
    //and writing data to file
    std::vector<double> efield;
    efield.push_back(E_ODE.initial_value());
    for (int i = 1; i < ITERATIONS; ++i) {
        double buffer = F_method_efield.step();
        solution << i * DX << " " << buffer << endl;
        efield.push_back(buffer);
    }
    solution.close();

    //instantiate potential field ODE and Method to solve it using the result
    //from the electric field analysis.
    PNJUNCTION::PotentialODE P_ODE(X0, H, efield, DX);
    PNJUNCTION::FORKMethod F_method_potential(DX, &P_ODE, X0);

    //creating output file and writing data to it.
    ofstream solutionP("solution_potential.txt");
    solutionP << X0 << " " << P_ODE.initial_value() << endl;
    for (int i = 1; i < ITERATIONS; ++i) {
        solutionP << i * DX << " " << F_method_potential.step() << endl;
    }
    solutionP.close();

    //creating exact solutions for comparison
    ofstream solutionPE("solution_potential_exact.txt");
    for (int i = 0; i < ITERATIONS; ++i) {
        solutionPE << i * DX << " " << P_ODE.exact_solution(i * DX) << endl;
    }
    solutionPE.close();
    ofstream solutionEE("solution_electric_exact.txt");
    for (int i = 0; i < ITERATIONS; ++i) {
        solutionEE << i * DX << " " << E_ODE.exact_solution(i * DX) << endl;
    }
    solutionEE.close();


    /* ****     Fourier Transform analysis      **** */

    // N bins
    int N=500 ;
    fftw_complex *rho, *rho_ft, *e_ft, *e;
    fftw_plan p;

    //creating arrays and output files to hold solutions.
    ofstream fft("solution_electric_fft.txt");
    rho = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    rho_ft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    e = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    e_ft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    //plan for forward fourier transforming the charge density
    p = fftw_plan_dft_1d(N, rho, rho_ft, FFTW_FORWARD, FFTW_ESTIMATE);

    for(int n=0; n < N ; ++n) { 
        rho[n][0] = E_ODE.first_derivative(0.0, n * XMAX/(double)N);
        rho[n][1] = 0.0;
    }

    //Fourier transforming the charge density
    fftw_execute(p);  

    // divide rho_ft by ik_n
    for(int n=0; n < N ; ++n) { 
        e_ft[n][0] = rho_ft[n][1]/(PI * n);
        e_ft[n][1] = - rho_ft[n][0] / (PI * n);
    }
    e_ft[0][0] = 0;
    e_ft[0][1] = 0;

    // Inverse Fourier transform of the Electric Field
    p = fftw_plan_dft_1d(N, e_ft, e, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    // writing the normalised values to file
    for(int n=0; n < N ; ++n) { 
        // normalisation
        double norm = H/(((e[N/2][0]/(double)N) - (e[0][0]/(double)N)));
        // y-shift
        double shift = e[0][0]/(double)N;

        fft << n * (XMAX/(double)N) << " " << norm * ((e[n][0]/(double)N) - shift) << endl;
    }
    fft.close(); 

    //Clean up
    fftw_destroy_plan(p);
    fftw_free(rho); fftw_free(rho_ft);

    return 0;

}
