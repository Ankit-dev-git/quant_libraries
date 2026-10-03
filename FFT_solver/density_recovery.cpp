#include "density_recovery.h"
#include <algorithm>

//----------------------------------------------------
// Recover Density
//----------------------------------------------------

std::vector<double> RecoverDensity(CF cf, const std::vector<double>& x, int N)
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