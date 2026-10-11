#include "grid_builder.h"
#include <stdexcept>

fourier::GridParams fourier::build_grid(const fourier::CharacteristicFuncResult& characteristic_func_XT, double K, double L, int N){
     if (N<=0 || K<=0 || L<=0 || characteristic_func_XT.sigma_terminal<=0)
        throw std::invalid_argument("N, L, K and sigma should be positive!");
    if (!std::isfinite(characteristic_func_XT.mu_terminal) || !std::isfinite(L) || !std::isfinite(K) || !std::isfinite(characteristic_func_XT.sigma_terminal))
        throw std::invalid_argument("L, K, mu and sigma should be finite!");       
    fourier::GridParams params;
    double dx = 2*L*characteristic_func_XT.sigma_terminal/N;
    params.du = 2*fourier::PI/(N * dx);
    int y = std::round((characteristic_func_XT.mu_terminal   - std::log(K))/dx);
    params.x_min = std::log(K) - (N/2 - y)*dx;
    params.points_per_sd = characteristic_func_XT.sigma_terminal/dx;
    params.n_dx_xmin = N/2 - y ;
    return params;
}