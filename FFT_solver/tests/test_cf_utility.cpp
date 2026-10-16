#include <gtest/gtest.h>
#include "characteristic_func.h"
#include <utility>
#include <ostream>
#include "FFT.h"

class CfGeneric: public testing::TestWithParam<fourier::CharacteristicFuncResult>{};
namespace{
    fourier::BlackScholes_params bs_params {.S_0 = 100, .r = 0.05, .sigma = 0.2, .T = 1};
    fourier::HestonParams heston_params {.S_0 = 100, .T = 1, .r = 0.0, .kappa = 1.5768, .gamma = 0.5751, .v_bar = 0.0398, .v0 = 0.0175, .rho = -0.5711};
}

namespace fourier{
    void PrintTo(const CharacteristicFuncResult& value, std::ostream* os){
        *os<< value.model;
    }
}

INSTANTIATE_TEST_SUITE_P(
    TestInputs,
    CfGeneric,
    testing::Values(fourier::charac_func_BS(bs_params), fourier::charac_func_Heston(heston_params)),
    [](const testing::TestParamInfo<fourier::CharacteristicFuncResult>& info){
        return info.param.model;
    }
);

TEST_P(CfGeneric, Phi0){
    fourier::CharacteristicFuncResult result = GetParam();
    Complex charac_func_0 = result.charac_func(0);
    double phi_0_real = charac_func_0.real();
    double phi_0_img = charac_func_0.imag();
    EXPECT_NEAR(phi_0_real, 1.0, 1e-14);
    EXPECT_NEAR(phi_0_img, 0.0, 1e-14);
}

TEST_P(CfGeneric, PhiUModulus){
    fourier::CharacteristicFuncResult result = GetParam();
    fourier::CF phi = result.charac_func;
    std::vector<double> u_vec {0.1, 1.0, 5.0};
    for (auto u: u_vec)
        EXPECT_LE(std::abs(phi(u)), 1.0)<< "u = " << u;
}

TEST_P(CfGeneric, PhiUConj){
    fourier::CharacteristicFuncResult result = GetParam();
    fourier::CF phi = result.charac_func;
    std::vector<double> u_vec {0.1, 1.0, 5.0};
    for (auto u: u_vec){
        Complex phi_minus_u = phi(-u);
        Complex phi_u_conj = std::conj(phi(u));
        EXPECT_NEAR(phi_u_conj.real(), phi_minus_u.real(), 1e-14) << "u = " << u;
        EXPECT_NEAR(phi_u_conj.imag(), phi_minus_u.imag(), 1e-14) << "u = " << u;
    }
}

TEST_P(CfGeneric, DecayLargeU){
    fourier::CharacteristicFuncResult result = GetParam();
    fourier::CF phi = result.charac_func;
    std::vector<double> u_vec {10, 50, 100};
    for (std::size_t i=1; i<u_vec.size(); ++i)
        EXPECT_GT(std::abs(phi(u_vec[i-1])),std::abs(phi(u_vec[i]))) << "u = " << u_vec[i] << " u(-1) = " << u_vec[i-1];
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