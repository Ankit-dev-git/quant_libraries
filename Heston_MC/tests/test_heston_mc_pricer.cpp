#include <gtest/gtest.h>
#include "qe_scheme.h"
#include "qe_pricer.h"
#include "modelparam_structs.h"
#include "normal.h"
#include "fourier_pricer.h"
#include "characteristic_func.h"
#include "grid_builder.h"
#include "BS_price.h"
#include <stdexcept>
#include <limits>
#include <vector>
#include <numeric>
#include <cmath>

TEST(HestonMCPricer, FourierCrossCheck){
// Robust case as the parameters violate Feller condition - 2κv̄ ≈ 0.16 < γ² ≈ 0.33 
    model_params::HestonParams heston_params {.S_0 = 100, .T = 1, .r = 0.05, .kappa = 2, .gamma = 0.5751, .v_bar = 0.0398, .v0 = 0.0398, .rho = -0.7};
    int N = 4096;
    double K = 100.0;
    fourier::CharacteristicFuncResult result_heston = fourier::charac_func_Heston(heston_params);
    fourier::GridParams grid_param = fourier::build_grid(result_heston, K, fourier::default_L, N);
    fourier::DensityResult heston_density = fourier::RecoverDensity(result_heston.charac_func, grid_param.du, grid_param.x_min, N);
    double put_price = fourier::put_price_from_density(heston_density, K, heston_params.r, heston_params.T);
    double call_price = put_price + heston_params.S_0 - K*std::exp(-heston_params.r*heston_params.T);

    heston_mc::HestonMCPricerResult result {heston_mc::heston_mc_pricer(heston_params, K, 100000, 50)};

    EXPECT_GT(result.std_error, 0.0) << "Heston MC std error: " << result.std_error;
    EXPECT_NEAR(call_price, result.price, 3*result.std_error) << "Fourier price: " << call_price << ", Heston MC price: " << result.price;
}

TEST(HestonMCPricer, MartingaleCheckDiscountedS){
// Robust case as the parameters violate Feller condition - 2κv̄ ≈ 0.16 < γ² ≈ 0.33 
    model_params::HestonParams heston_params {.S_0 = 100, .T = 1, .r = 0.05, .kappa = 2, .gamma = 0.5751, .v_bar = 0.0398, .v0 = 0.0398, .rho = -0.7};

// At K = 0, price is E[e^xT] * e^-rt, which is discounted value of Terminal price expectation
    heston_mc::HestonMCPricerResult result = heston_mc::heston_mc_pricer(heston_params, 0.0, 100000);
    double discounted_expected_sT = result.price;
    EXPECT_GT(result.std_error, 0.0);
    EXPECT_NEAR(discounted_expected_sT, heston_params.S_0, 3*result.std_error);
}

TEST(HestonMCPricer, BSReduction){
// At gamma = 0.01, heston impact would come from rho*gamma + O(gamma^2), which is smaller compared to 3*S.E (0.14). 
// As v0=vbar and gamma is small, the kappa factors impact is zero and the model reduced to BS.
    model_params::HestonParams heston_params {.S_0 = 100, .T = 10, .r = 0.05, .kappa = 2, .gamma = 0.01, .v_bar = 0.0398, .v0 = 0.0398, .rho = -0.7};
    double K1 = 100;
    heston_mc::HestonMCPricerResult result1 = heston_mc::heston_mc_pricer(heston_params, K1, 100000);
// The 0 is for initial time t0 which is zero in our case; 
    double bs_price1 = bs::option_price(OptionType::call, heston_params.S_0, K1, std::sqrt(heston_params.v0), 0, heston_params.T, heston_params.r);
    EXPECT_GT(result1.std_error, 0.0);
    EXPECT_NEAR(result1.price, bs_price1, result1.std_error*3);
}

