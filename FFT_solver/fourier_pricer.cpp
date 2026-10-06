#include "fourier_pricer.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

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