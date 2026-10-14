#include <gtest/gtest.h>
#include "characteristic_func.h"
#include <utility>

TEST(CfGeneric, Phi0){
    fourier::BlackScholes_params param {.S_0 = 100, .r = 0.05, .sigma = 0.2, .T = 1};
    fourier::CharacteristicFuncResult result1 = fourier::charac_func_BS(param);
    double phi_0_real = result1.charac_func(0).real();
    double phi_0_img = result1.charac_func(0).imag();
    EXPECT_NEAR(phi_0_real, 1.0, 1e-14);
    EXPECT_NEAR(phi_0_img, 0.0, 1e-14);
}
TEST(CfGeneric, PhiUModulus){
    fourier::BlackScholes_params param {.S_0 = 100, .r = 0.05, .sigma = 0.2, .T = 1};
    fourier::CharacteristicFuncResult result1 = fourier::charac_func_BS(param);
    fourier::CF phi = result1.charac_func;
    std::vector<double> u_vec {0.1, 1.0, 5.0};
    for (auto u: u_vec)
        EXPECT_LE(std::abs(phi(u)), 1.0)<< "u = " << u;
}
TEST(CfGeneric, PhiUConj){
    fourier::BlackScholes_params param {.S_0 = 100, .r = 0.05, .sigma = 0.2, .T = 1};
    fourier::CharacteristicFuncResult result1 = fourier::charac_func_BS(param);
    fourier::CF phi = result1.charac_func;
    std::vector<double> u_vec {0.1, 1.0, 5.0};
    for (auto u: u_vec){
        Complex phi_minus_u = phi(-u);
        EXPECT_NEAR(phi(u).real(), std::conj(phi_minus_u).real(), 1e-14) << "u = " << u;
        EXPECT_NEAR(phi(u).imag(), std::conj(phi_minus_u).imag(), 1e-14) << "u = " << u;
    }
}

TEST(CfGeneric, CfCumulants){
    fourier::BlackScholes_params param {.S_0 = 100, .r = 0.05, .sigma = 0.2, .T = 2};
    fourier::CharacteristicFuncResult result = fourier::charac_func_BS(param);
    double mu_T = std::log(param.S_0) + (param.r - param.sigma*param.sigma/2)*param.T;
    double var_T = param.sigma*param.sigma*param.T;
    std::pair<double, double> mean_var {fourier::charac_func_cumulants(result.charac_func)};
    EXPECT_NEAR(mu_T, mean_var.first, 1e-12);
// BS log(phi) has no truncation error, thus this check doesn't check that. Heston reducing to BS and then calculating cumulants will check that;
    EXPECT_NEAR(var_T, mean_var.second, 1e-9);
}

TEST(CfGeneric, PosNegH){
    fourier::BlackScholes_params param {.S_0 = 100, .r = 0.05, .sigma = 0.2, .T = 2};
    fourier::CharacteristicFuncResult result = fourier::charac_func_BS(param);
    std::pair<double, double> mean_var_pos_h {fourier::charac_func_cumulants(result.charac_func, 1e-3)};
    std::pair<double, double> mean_var_neg_h {fourier::charac_func_cumulants(result.charac_func, -1e-3)};
    EXPECT_NEAR(mean_var_pos_h.first, mean_var_neg_h.first, 1e-12);
    EXPECT_NEAR(mean_var_pos_h.second, mean_var_neg_h.second, 1e-9);
}