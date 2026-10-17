#include "grid_builder.h"
#include <stdexcept>

// This function takes as input the CharacteristicFuncResult struct object, from which it extract c1, c2 and c4; It also takes K, L and N as user input;
// This function returns the Grid_parameters, mainly du and x_min to ensure lnK lies on the grid and there is no wrap-around. 
// This parameters are then passed directly into RecoverDensity function. points_per_sd and n_dx_xmin are not used in pricing and only used for testing purposes.
fourier::GridParams fourier::build_grid(const fourier::CharacteristicFuncResult& characteristic_func_XT, double K, double L, int N){
     if (N<=0 || K<=0 || L<=0 || characteristic_func_XT.sigma_terminal<=0)
        throw std::invalid_argument("N, L, K and sigma should be positive!");
    if (!std::isfinite(characteristic_func_XT.mu_terminal) || !std::isfinite(L) || !std::isfinite(K) || !std::isfinite(characteristic_func_XT.sigma_terminal) || !std::isfinite(characteristic_func_XT.c4))
        throw std::invalid_argument("L, K, mu, sigma and c4 should be finite!");       

// Adjusted for kurtosis by adjsuting the effective sigma
    double c4 = (characteristic_func_XT.c4<0) ? 0:characteristic_func_XT.c4;
    double s_eff = (std::sqrt(characteristic_func_XT.sigma_terminal*characteristic_func_XT.sigma_terminal + std::sqrt(c4)));
    fourier::GridParams params;
    double dx = 2*L*s_eff/N;
    params.du = 2*fourier::PI/(N * dx);
    int y = std::lround((characteristic_func_XT.mu_terminal   - std::log(K))/dx);
    params.x_min = std::log(K) - (N/2 - y)*dx;
    params.points_per_sd = characteristic_func_XT.sigma_terminal/dx;
    params.n_dx_xmin = N/2 - y ;
    return params;
}