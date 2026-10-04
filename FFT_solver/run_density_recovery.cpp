#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <numbers>
#include "FFT.h"
#include "normal.h"
#include "density_recovery.h"

int main()
{
    double mu = 0.0;
    double sigma = 1.0;

    int N = 64; // same for both X and u
    double u_span = 20.0; 
    double du = u_span/N;
    /* grid is from [0, 20-delta_u]; From 20, the grid repeats
    u_min = 0; u_span = (N)*delta_u, thus delta_u = 20/64; 
    We know delta_x*delta_u = 2*pi/N = pi/32 (constant as N is fixed);  
    delta_x = pi/(32*delta_u) = (pi/32)*(64/20); */


// Characteristic Function of Normal Distribution
    CF normal_cf = [=](double u)
        {
            Complex i(0.0, 1.0);

            return std::exp(i * mu * u - 0.5 * sigma * sigma * u * u);
        };


// Build calculating x_min and dx
    double dx = 2.0 * PI / (du*N);
    double x_min = -N * dx / 2.0; 

    /* x_min = −(N/2)·dx; so X runs from −N·dx/2 to (N·dx/2 − dx); 
    It's centred, but asymmetric by one step: 
    for N = 64, from −10.05 to +9.74. */


// Calculating density using fourier inversion and DFT.
    DensityResult recovered_density_obj = RecoverDensity(normal_cf, du, x_min, N);

//------------------------------------------------
// Compute exact pdf for testing and comparison
//------------------------------------------------
    std::vector<double> exact_pdf(N);
    
    for (int k = 0; k < N; ++k)
    {
        exact_pdf[k] = normal::norm_pdf(recovered_density_obj.x[k], mu, sigma);
    }

    //------------------------------------------------
    // Export data for plotting
    //------------------------------------------------

    std::ofstream file("density_comparison.csv");

    file << "x,recovered_pdf,exact_pdf\n";

    for (int k = 0; k < N; ++k)
    {
        file
            << recovered_density_obj.x[k] << ","
            << recovered_density_obj.density[k] << ","
            << exact_pdf[k] << "\n";
    }

    file.close();

    std::cout
        << "Data exported to density_comparison.csv"
        << std::endl;

    return 0;
}