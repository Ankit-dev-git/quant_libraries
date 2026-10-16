#include <gtest/gtest.h>
#include "characteristic_func.h"
#include "FFT.h"
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
    EXPECT_TRUE(std::isfinite(mean));
    EXPECT_TRUE(std::isfinite(sigma));
    EXPECT_GE(sigma, 0.0 + 1e-12);
    EXPECT_NEAR(sigma*sigma, fourier::charac_func_cumulants(result.charac_func).second, 1e-12);
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