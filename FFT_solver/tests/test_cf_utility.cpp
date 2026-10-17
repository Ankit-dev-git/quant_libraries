#include <gtest/gtest.h>
#include "characteristic_func.h"
#include <utility>
#include <ostream>
#include <algorithm>
#include "FFT.h"
#include "grid_builder.h"
#include "density_recovery.h"
#include "fourier_pricer.h"

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

TEST_P(CfGeneric, StartEndDensityValueZero){
    int N = 1024;
    double K = 100.0;
    fourier::CharacteristicFuncResult result = GetParam();
    fourier::GridParams grid_param = fourier::build_grid(result, K, fourier::default_L, N);
    std::vector<double> density = fourier::RecoverDensity(result.charac_func, grid_param.du, grid_param.x_min, N).density;
    double peak = (*std::max_element(density.begin(), density.end()));
    EXPECT_NEAR(density[0], 0.0, 1e-6 * peak);
    EXPECT_NEAR(density.back(), 0.0, 1e-6* peak);
}

TEST(CfGeneric, CfCumulants){
    fourier::BlackScholes_params param {.S_0 = 100, .r = 0.05, .sigma = 0.2, .T = 2};
    fourier::CharacteristicFuncResult result = fourier::charac_func_BS(param);
    double mu_T = std::log(param.S_0) + (param.r - param.sigma*param.sigma/2)*param.T;
    double var_T = param.sigma*param.sigma*param.T;
    fourier::Cumulants cumulants {fourier::charac_func_cumulants(result.charac_func)};
    EXPECT_NEAR(mu_T, cumulants.c1, 1e-12);
// BS log(phi) has no truncation error, thus this check doesn't check that. Heston reducing to BS and then calculating cumulants will check that;
    EXPECT_NEAR(var_T, cumulants.c2, 1e-8);
    EXPECT_NEAR(0.0, cumulants.c4, 1e-6);
}

TEST(CfGeneric, CfCumulantsSyntheticC4){
    // log φ(u) is an exact cumulant polynomial, so the true c1, c2, c4 are known.
    // Not a valid distribution (|φ| grows at large u), but the estimator only evaluates u ≤ 2e-2.
    const double c1 = 4.6, c2 = 0.08, c3 = -0.01, c4 = 0.005, c6 = 0.001;
    fourier::CF synthetic_cf = [=](double u){
        Complex i(0.0, 1.0);
        Complex log_phi = i*c1*u - c2*u*u/2.0 - i*c3*u*u*u/6.0
                        + c4*u*u*u*u/24.0 - c6*u*u*u*u*u*u/720.0;
        return std::exp(log_phi);
    };
    fourier::Cumulants cumulants = fourier::charac_func_cumulants(synthetic_cf);
    // c1: truncation c3·h²/6 ~ 2e-11 at h = 1e-4
    EXPECT_NEAR(cumulants.c1, c1, 1e-10);
    // c2: round-off ~2ε/h² ~ 2e-8; truncation c4·h²/12 ~ 4e-12
    EXPECT_NEAR(cumulants.c2, c2, 1e-7);
    // c4: round-off ~16ε/h⁴ ~ 1e-7; truncation c6·h²/6 ~ 2e-8 at h = 1e-2
    EXPECT_NEAR(cumulants.c4, c4, 1e-6);
}

TEST(CfGeneric, DuWidening){
    fourier::BlackScholes_params param {.S_0 = 100, .r = 0.05, .sigma = 0.2, .T = 2};
    fourier::CharacteristicFuncResult result = fourier::charac_func_BS(param);
    result.c4 = 0.2;
    fourier::GridParams grid_param =  fourier::build_grid(result, 100);
//  du = π/(L·√(c₂ + √c₄))
    double du_calc = fourier::PI/(fourier::default_L*(std::sqrt(result.sigma_terminal*result.sigma_terminal + std::sqrt(result.c4))));
    EXPECT_NEAR(grid_param.du, du_calc , 1e-12);
}

TEST(CfGeneric, PosNegH){
    fourier::BlackScholes_params param {.S_0 = 100, .r = 0.05, .sigma = 0.2, .T = 2};
    fourier::CharacteristicFuncResult result = fourier::charac_func_BS(param);
    fourier::Cumulants cumulants_pos_h {fourier::charac_func_cumulants(result.charac_func, 1e-3)};
    fourier::Cumulants cumulants_neg_h {fourier::charac_func_cumulants(result.charac_func, -1e-3)};
    EXPECT_NEAR(cumulants_pos_h.c1, cumulants_neg_h.c1, 1e-12);
    EXPECT_NEAR(cumulants_pos_h.c2, cumulants_neg_h.c2, 1e-9);
}