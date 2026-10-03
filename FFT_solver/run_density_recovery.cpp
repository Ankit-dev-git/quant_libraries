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

    int N = 64;
    double u_max = 20.0;

    //------------------------------------------------
    // Characteristic Function of Normal Distribution
    //------------------------------------------------

    CF normal_cf = [=](double u)
        {
            Complex i(0.0, 1.0);

            return std::exp(i * mu * u - 0.5 * sigma * sigma * u * u);
        };


    //------------------------------------------------
    // Build x-grid
    //------------------------------------------------

    double dx = 2.0 * PI / u_max;

    std::vector<double> x_grid(N);

    for (int k = 0; k < N; ++k)
        x_grid[k] = -N * dx / 2.0 + k * dx;

    //------------------------------------------------
    // Compute exact PDF
    //------------------------------------------------

    std::vector<double> recovered_density = RecoverDensity(normal_cf, x_grid, N);

    std::vector<double> exact_pdf(N);

    for (int k = 0; k < N; ++k)
    {
        exact_pdf[k] = normal::norm_pdf(x_grid[k], mu, sigma);
    }

    //------------------------------------------------
    // Export data for plotting
    //------------------------------------------------

    std::ofstream file("density_comparison.csv");

    file << "x,recovered_pdf,exact_pdf\n";

    for (int k = 0; k < N; ++k)
    {
        file
            << x_grid[k] << ","
            << recovered_density[k] << ","
            << exact_pdf[k] << "\n";
    }

    file.close();

    std::cout
        << "Data exported to density_comparison.csv"
        << std::endl;

    return 0;
}