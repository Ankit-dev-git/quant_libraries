#include "fourier_pricer.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

// This function calculates call price by taking just characteristic function and x-grid (passed through result struct) and K,r and T;
// Direct call integration: sum of max(e^x - K, 0) * f(x) * dx over the full grid.
// Accurate when the window's right edge is moderate (e.g. T = 1 reference case).
// For wide windows (long T / fat tails, large s_eff), e^x amplifies the density's
// error floor at the right edge (measured: E[S_T] off by ~3.4e5 at T = 10, x_max = 27).
// In that regime use put_price_from_density + put-call parity:
// C = P + S0 - K * exp(-r*T)   (put payoff is bounded by K, so immune to the right tail).
double fourier::call_price_from_density(const DensityResult& result, double K, double r, double T){
    double price {0};
    if (K<0 || result.x.size()<2)
        throw std::invalid_argument("Strike Price K should not be less than zero and grid length greater than 1");
    
    double x_width = result.x[1] - result.x[0];

//  Calculating the call price. Note that density contains pdf and not probabilities, so we need to multiply with dx  for sum to be area.
    for (size_t i =0; i<result.x.size(); ++i)
        price+= std::max(std::exp(result.x[i])-K,0.0)*result.density[i]*x_width;
    return price*std::exp(-1*r*T);
}

double fourier::put_price_from_density(const DensityResult& result, double K, double r, double T){
    double price {0};
    if (K<0 || result.x.size()<2)
        throw std::invalid_argument("Strike Price K should not be less than zero and grid length greater than 1");
    
    double x_width = result.x[1] - result.x[0];

//  Calculating the call price. Note that density contains pdf and not probabilities, so we need to multiply with dx  for sum to be area.
    for (size_t i =0; i<result.x.size(); ++i)
        price+= std::max(K-std::exp(result.x[i]),0.0)*result.density[i]*x_width;
//  Terminal price is calculated from above; We return current price by multiplying with discount factor;
    return price*std::exp(-r*T);
}