#include "FFT.h"
#include "normal.h"
#include "characteristic_func.h"
#include "density_recovery.h"
#include "fourier_pricer.h"
#include "BS_price.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <numbers>
#include <stdexcept>

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
    fourier::CF normal_cf = [=](double u)
        {
            Complex i(0.0, 1.0);

            return std::exp(i * mu * u - 0.5 * sigma * sigma * u * u);
        };


// Build calculating x_min and dx
    double dx = 2.0 * fourier::PI / (du*N);
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

// Check for BS_charac_func
    fourier::BlackScholes_params params {.S_0 = 100, .r = 0.05, .sigma = 0.2, .T = 1};
    fourier::BSInvFourResult BS_cf_mu = fourier::charac_func_BS(params);
    std::cout<<"Black_scholes characteristic function at u=0: "<< BS_cf_mu.charac_func(0) << std::endl;
    std::cout<<"Black_scholes characteristic function absolute value at u=10 and T= "<< params.T <<": "<< std::abs(BS_cf_mu.charac_func(10)) << std::endl;

    N = 64; // same for both X and u
    u_span = 100.0;
    double K = 100.0;
    double analytical_price = bs::option_price(OptionType::call, params.S_0, K, params.sigma, 0.0, params.T, params.r);
    double last_price {0.0};
    std::cout << "Analytical BS Call Price: " <<  analytical_price << ";" << std::endl; 
// Checking if Inverse fourier FFT works
    for (int i=1; i<=10; i++){
    // Error decrease by 4 for doubling of dx.
        du = u_span/N;
        dx = 2.0 * fourier::PI / (du*N);
        int dx_shift = static_cast<int>((BS_cf_mu.mu - std::log(K))/dx);
        double x_min_bs = std::log(K) - (N/2 - dx_shift)*dx;
        // double x_min_bs = BS_cf_mu.mu - N/2*dx;
        DensityResult bs_recovered_density_obj = RecoverDensity(BS_cf_mu.charac_func, du, x_min_bs, N);
        double fourier_call_price {fourier::call_price_from_density(bs_recovered_density_obj, K, params.r, params.T)};
        
        // std::cout<< "Fourier BS call price: "<< fourier_call_price << std::endl; 
        std::cout << "N: "<< N<< "; Direct error: " << analytical_price - fourier_call_price << std::endl;
    
        // Richardson Method to make error decrease by 2^4 for every doubling of dx
        if (last_price){
            double method_2_price = (4*fourier_call_price - last_price)/3;
            std::cout << "N: "<< N<< "; Richardson error: " << analytical_price - method_2_price << std::endl;
        }
        last_price = fourier_call_price;
        u_span*=2;
        N*=2;
    }

    return 0;
} 