#include <gtest/gtest.h>
#include "characteristic_func.h"
#include "FFT.h"
#include "grid_builder.h"
#include "density_recovery.h"
#include "fourier_pricer.h"
#include <stdexcept>
#include <limits>
#include <vector>

TEST(HestonCf, ThrowsOnBadInputs){
    fourier::HestonParams no_throw_params {.S_0 = 100, .T = 1, .r = 0.0, .kappa = 1.5768, .gamma = 0.5751, .v_bar = 0.0398, .v0 = 0.0175, .rho = -0.5711};
    EXPECT_NO_THROW(fourier::charac_func_Heston(no_throw_params));

    std::vector<double> vec {0.0, -5.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    
    fourier::HestonParams s0_throw_params = no_throw_params;
    for (double& elem:vec){
        s0_throw_params.S_0 = elem;
        EXPECT_THROW(fourier::charac_func_Heston(s0_throw_params), std::invalid_argument)<< "elem = " << elem;
    }

    fourier::HestonParams t_throw_params = no_throw_params;
    for (double& elem:vec){
        t_throw_params.T = elem;
        EXPECT_THROW(fourier::charac_func_Heston(t_throw_params), std::invalid_argument)<< "elem = " << elem;
    }

    fourier::HestonParams kappa_throw_params = no_throw_params;
    for (double& elem:vec){
        kappa_throw_params.kappa = elem;
        EXPECT_THROW(fourier::charac_func_Heston(kappa_throw_params), std::invalid_argument)<< "elem = " << elem;
    }

    fourier::HestonParams gamma_throw_params = no_throw_params;
    for (double& elem:vec){
        gamma_throw_params.gamma = elem;
        EXPECT_THROW(fourier::charac_func_Heston(gamma_throw_params), std::invalid_argument)<< "elem = " << elem;
    }

    fourier::HestonParams vbar_throw_params = no_throw_params;
    for (double& elem:vec){
        vbar_throw_params.v_bar = elem;
        EXPECT_THROW(fourier::charac_func_Heston(vbar_throw_params), std::invalid_argument)<< "elem = " << elem;
    }

    std::vector<double> vec_r {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    fourier::HestonParams r_throw_params = no_throw_params;
    for (double& elem:vec_r){
        r_throw_params.r = elem;
        EXPECT_THROW(fourier::charac_func_Heston(r_throw_params), std::invalid_argument)<< "elem = " << elem;
    }

    std::vector<double> vec_v0 {-5, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    fourier::HestonParams v0_throw_params = no_throw_params;
    for (double& elem:vec_v0){
        v0_throw_params.v0 = elem;
        EXPECT_THROW(fourier::charac_func_Heston(v0_throw_params), std::invalid_argument)<< "elem = " << elem;
    }

    std::vector<double> vec_rho {-5, +2, -1.0, +1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    fourier::HestonParams rho_throw_params = no_throw_params;
    for (double& elem:vec_rho){
        rho_throw_params.rho = elem;
        EXPECT_THROW(fourier::charac_func_Heston(rho_throw_params), std::invalid_argument)<< "elem = " << elem;
    }
}

TEST(HestonCf, AcceptsZeroV0){
    fourier::HestonParams no_throw_params {.S_0 = 100, .T = 1, .r = 0.0, .kappa = 1.5768, .gamma = 0.5751, .v_bar = 0.0398, .v0 = 0.0, .rho = -0.5711};
    EXPECT_NO_THROW(fourier::charac_func_Heston(no_throw_params));
}

TEST(HestonCf, AcceptsNegR){
    fourier::HestonParams no_throw_params {.S_0 = 100, .T = 1, .r = -0.5, .kappa = 1.5768, .gamma = 0.5751, .v_bar = 0.0398, .v0 = 0.03, .rho = -0.5711};
    EXPECT_NO_THROW(fourier::charac_func_Heston(no_throw_params));
}

TEST(HestonCf, ModelName){
    fourier::HestonParams no_throw_params {.S_0 = 100, .T = 1, .r = 0.0, .kappa = 1.5768, .gamma = 0.5751, .v_bar = 0.0398, .v0 = 0.0175, .rho = -0.5711};
    EXPECT_EQ(fourier::charac_func_Heston(no_throw_params).model, "Heston");
}

TEST(HestonCf, TerminalMoments){
    fourier::HestonParams no_throw_params {.S_0 = 100, .T = 1, .r = 0.0, .kappa = 1.5768, .gamma = 0.5751, .v_bar = 0.0398, .v0 = 0.0175, .rho = -0.5711};
    fourier::CharacteristicFuncResult result = fourier::charac_func_Heston(no_throw_params);
    double mean = result.mu_terminal;
    double sigma = result.sigma_terminal;
    double c4 = result.c4;
    EXPECT_TRUE(std::isfinite(mean));
    EXPECT_TRUE(std::isfinite(sigma));
    EXPECT_TRUE(std::isfinite(c4));
    EXPECT_GE(sigma, 0.0 + 1e-12);
    EXPECT_GE(c4, 0.0 + 1e-12);
    EXPECT_NEAR(sigma*sigma, fourier::charac_func_cumulants(result.charac_func).c2, 1e-12);
}

TEST(HestonCf, DecayLargeU){
    fourier::HestonParams heston_params {.S_0 = 100, .T = 1, .r = 0.0, .kappa = 1.5768, .gamma = 0.5751, .v_bar = 0.0398, .v0 = 0.0175, .rho = -0.5711};
    fourier::CharacteristicFuncResult result = fourier::charac_func_Heston(heston_params);
    fourier::CF phi = result.charac_func;
// Heston φ decays exponentially in u (≈ e^{−0.114u} for these params), not Gaussian like BS, hence the looser threshold / larger u.
    std::vector<double> u_vec {150, 200, 300};
    for (auto u: u_vec)
        EXPECT_LT(std::abs(phi(u)), 1e-5)<< "u = " << u;
}

TEST(HestonCf, ReducedToBS){
    fourier::HestonParams heston_params {.S_0 = 100, .T = 1, .r = 0.0, .kappa = 1.5768, .gamma = 0.01, .v_bar = 0.04, .v0 = 0.04, .rho = -0.5711};
    fourier::CharacteristicFuncResult result_heston = fourier::charac_func_Heston(heston_params);
    fourier::CF phi_heston = result_heston.charac_func;
    fourier::BlackScholes_params bs_params {.S_0 = 100, .r = 0.0, .sigma = 0.2, .T = 1};
    fourier::CharacteristicFuncResult result_bs = fourier::charac_func_BS(bs_params);
    fourier::CF phi_bs = result_bs.charac_func;
    std::vector<double> u_vec {0.5, 1, 2, 5};
    EXPECT_NEAR(fourier::charac_func_cumulants(phi_heston).c1, fourier::charac_func_cumulants(phi_bs).c1, 1e-3);
    EXPECT_NEAR(fourier::charac_func_cumulants(phi_heston).c2, fourier::charac_func_cumulants(phi_bs).c2, 1e-3);
    EXPECT_NEAR(fourier::charac_func_cumulants(phi_heston).c4, fourier::charac_func_cumulants(phi_bs).c4, 1e-3);
// measured max error 2.5e-3 at u = 5, γ = 0.01; grows with u, shrinks with γ
    for (const auto& u: u_vec){
        EXPECT_NEAR(phi_heston(u).imag(), phi_bs(u).imag(), 5e-3)<< "u = "<< u;
        EXPECT_NEAR(phi_heston(u).real(), phi_bs(u).real(), 5e-3)<< "u = "<< u;
    }
    std::vector<double> gamma_vec = {1e-3, 1e-4};
    for (const auto& gamma: gamma_vec){
        heston_params.gamma = gamma;
        fourier::CharacteristicFuncResult result_heston = fourier::charac_func_Heston(heston_params);
        fourier::CF phi_heston_upd = result_heston.charac_func;
        for (const auto& u: u_vec){
            double error_1 = std::abs(phi_bs(u) - phi_heston(u));
            double error_2 = std::abs(phi_bs(u) - phi_heston_upd(u));
            EXPECT_GE(error_1/5, error_2);
        }
        phi_heston = phi_heston_upd;
    }
}

TEST(HestonCf, ContinuityLongMaturity){
    fourier::HestonParams heston_params {.S_0 = 100, .T = 10, .r = 0.0, .kappa = 1, .gamma = 1, .v_bar = 0.09, .v0 = 0.09, .rho = -0.9};
    fourier::CharacteristicFuncResult result_heston = fourier::charac_func_Heston(heston_params);
    fourier::CF phi_heston = result_heston.charac_func;
    double u {0}, mean_XT = fourier::charac_func_cumulants(phi_heston).c1, diff;
    Complex phi_u_trap_part, i(0.0, 1), phi_um1_trap_part {phi_heston(0)*std::exp(-i*u*mean_XT)};
    
    for (int k =1; k<=1000; k++){
        u = k*0.05;
        phi_u_trap_part = phi_heston(u)* std::exp(-i*u*mean_XT);
        diff = std::abs(phi_u_trap_part - phi_um1_trap_part);
        EXPECT_LE(diff, 0.05);
        phi_um1_trap_part = phi_u_trap_part;
    }
}

TEST(HestonCf, FellerViolatedStillValid){
    fourier::HestonParams heston_params {.S_0 = 100, .T = 10, .r = 0.0, .kappa = 1, .gamma = 1, .v_bar = 0.09, .v0 = 0.09, .rho = -0.9};
//  2κv̄ = 0.18 < γ² = 1; Thus failing feller condition, but characteristic function method should be unaffected and throw no error;
    EXPECT_NO_THROW(fourier::charac_func_Heston(heston_params));
    fourier::CharacteristicFuncResult result_heston = fourier::charac_func_Heston(heston_params);
    fourier::CF phi_heston = result_heston.charac_func;
    std::vector<double> u_vec {0.5, 1, 2, 5};
    for (const auto& u: u_vec){
        EXPECT_TRUE(std::isfinite(std::abs(phi_heston(u))));
        EXPECT_LE(std::abs(phi_heston(u)), 1.0);
    }
    fourier::Cumulants cumulants = fourier::charac_func_cumulants(phi_heston);
    double mean_XT = cumulants.c1;
    double var_XT = cumulants.c2;
    EXPECT_TRUE(std::isfinite(mean_XT));
    EXPECT_TRUE(std::isfinite(var_XT));
    EXPECT_GT(var_XT, 0.0);
}

TEST(HestonCf, ReferencePrice){
    int N = 1024;
    double K = 100.0;
    fourier::HestonParams heston_params {.S_0 = 100, .T = 1.0, .r = 0.0, .kappa = 1.5768, .gamma = 0.5751, .v_bar = 0.0398, .v0 = 0.0175, .rho = -0.5711};
    fourier::CharacteristicFuncResult result_heston = fourier::charac_func_Heston(heston_params);
    fourier::GridParams grid_param = fourier::build_grid(result_heston, K, fourier::default_L, N);
    fourier::DensityResult heston_density = fourier::RecoverDensity(result_heston.charac_func, grid_param.du, grid_param.x_min, N);
    double price = fourier::call_price_from_density(heston_density, K, heston_params.r, heston_params.T);
    EXPECT_NEAR(price, 5.785155, 1e-3);
}

TEST(HestonCf, LongMaturityT){
    int N = 4096;
    double K = 100;
    fourier::HestonParams heston_params {.S_0 = 100, .T = 10.0, .r = 0.0, .kappa = 1.0, .gamma = 1.0, .v_bar = 0.09, .v0 = 0.09, .rho = -0.9};
    fourier::CharacteristicFuncResult result_heston = fourier::charac_func_Heston(heston_params);
    fourier::GridParams grid_param = fourier::build_grid(result_heston, K, fourier::default_L, N);
    fourier::DensityResult heston_density = fourier::RecoverDensity(result_heston.charac_func, grid_param.du, grid_param.x_min, N);
    double put_price = fourier::put_price_from_density(heston_density, K, heston_params.r, heston_params.T);
//  1024 → 27.9013, 2048 → 27.9075, 4096 27.909; the wide T = 10 window (s_eff ≈ 2.9) needs more points.
    double call_price = put_price + heston_params.S_0 - K*std::exp(-heston_params.r*heston_params.T);
    EXPECT_NEAR(call_price, 27.91, 5e-3);
}
