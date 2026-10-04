#include "density_recovery.h"

//----------------------------------------------------
// Recover Density
//----------------------------------------------------

DensityResult RecoverDensity(CF cf, double du, double x_min, int N)
{
    Complex i(0.0, 1.0);

// Creating empty DensityResult object
    DensityResult result;

// Building the u-grid
    std::vector<double> u(N);
    for (int k = 0; k < N; ++k)
        u[k] = k * du; 
    /* u-grid: u_k = k*du, k = 0..N-1  ->  [0, (N-1)*du]
       The boundary correction below turns the sum into the trapezoid rule on [0, (N-1)*du];
       the integral beyond the last point is dropped (negligible when phi(u) has decayed). */

    result.x.resize(N);
// Building the x-grid
    double dx = 2.0 * PI / (du*N);
    for (int k = 0; k < N; ++k)
        result.x[k] = x_min + k * dx;


//------------------------------------------------
// Build the vector the DFT acts on (this is not characteristic function). 
//------------------------------------------------
    std::vector<Complex> phi(N);

//   phi[k] = e^{-i * x_min * u_k} * cf(u_k)
    for (int k = 0; k < N; ++k)
        phi[k] = std::exp(-i * x_min * u[k]) * cf(u[k]);

// Solving FFT using the phi vector created above
    FFT(phi); 
    // This function takes the phi by reference and alter it in place and return.
    // It multiplies the phi vector with M matrix (matrix of omega)
    // If we dont want to use FFT, simply do matrix multiplication of M and phi, however that takes O(N*N) instead of O(NlogN)

//------------------------------------------------
// Recover density: f(x_k) for k=0..,N-1
//------------------------------------------------
    result.density.resize(N);
    Complex gamma_1 = std::exp(-i * result.x[0] * u[0]) * cf(u[0]);
    Complex cf_last = cf(u[N - 1]);
    for (int k = 0; k < N; ++k){
    // Boundary correction
        Complex gamma_2 = std::exp(-i * result.x[k] * u[N - 1]) * cf_last;
        Complex phi_boundary = 0.5 * (gamma_1 + gamma_2);
        result.density[k] = du / PI * (phi[k] - phi_boundary).real();
    }

// Final Density object returned
    return result;
}