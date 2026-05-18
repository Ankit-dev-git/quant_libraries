#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <functional>
#include <fstream>
#include "fft.h"

using Complex = std::complex<double>;

const double PI = 3.14159265358979323846;

// CF signature
using CF = std::function<Complex(double)>;

//----------------------------------------------------
// Recover Density
//----------------------------------------------------

std::vector<double> RecoverDensity(CF cf, const std::vector<double>& x, int N = 64)
{
    Complex i(0.0, 1.0);
    double u_max = 20.0;
    double du = u_max / N;

    std::vector<double> u(N);

    for (int k = 0; k < N; ++k)
    {
        u[k] = k * du;
    }

// x-grid

    double b = *std::min_element(x.begin(), x.end());

    double dx = 2.0 * PI / (N * du);

    std::vector<double> x_i(N);

    for (int k = 0; k < N; ++k)
    {
        x_i[k] = b + k * dx;
    }

    //------------------------------------------------
    // Characteristic function samples
    //------------------------------------------------

    std::vector<Complex> phi(N);

    for (int k = 0; k < N; ++k)
    {
        phi[k] = std::exp(-i * b * u[k]) * cf(u[k]);
    }

    //------------------------------------------------
    // Boundary correction
    //------------------------------------------------

    Complex gamma_1 = std::exp(-i * x_i[0] * u[0]) * cf(u[0]);

    Complex gamma_2 = std::exp(-i * x_i[0] * u[N - 1]) * cf(u[N - 1]);

    Complex phi_boundary = 0.5 * (gamma_1 + gamma_2);

    FFT(phi);

    //------------------------------------------------
    // Recover density
    //------------------------------------------------

    std::vector<double> density(N);

    for (int k = 0; k < N; ++k)
    {
        density[k] = du / PI * (phi[k] - phi_boundary).real();
    }

    return density;
}

//----------------------------------------------------
// Exact Normal PDF
//----------------------------------------------------

double NormalPDF(double x, double mu, double sigma)
{
    double coeff = 1.0 / (sigma * std::sqrt(2.0 * PI));

    double exponent = -0.5 * std::pow((x - mu) / sigma, 2.0);

    return coeff * std::exp(exponent);
}


//----------------------------------------------------
// Main
//----------------------------------------------------

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
    {
        x_grid[k] =
            -N * dx / 2.0 + k * dx;
    }

    //------------------------------------------------
    // Compute exact PDF
    //------------------------------------------------

    std::vector<double> recovered_density = RecoverDensity(normal_cf, x_grid, N);

    std::vector<double> exact_pdf(N);

    for (int k = 0; k < N; ++k)
    {
        exact_pdf[k] = NormalPDF(x_grid[k], mu, sigma);
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