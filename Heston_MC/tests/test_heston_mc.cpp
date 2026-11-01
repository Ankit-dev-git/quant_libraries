#include <gtest/gtest.h>
#include "qe_scheme.h"
#include "normal.h"
#include <stdexcept>
#include <limits>
#include <vector>
#include <numeric>
#include <cmath>

TEST(HestonMC, MomentMatchQuardatic){
    int N {100000};
    std::vector<double> u_vec, vt_vec{};
    for (int i=1; i<=N; ++i)
        u_vec.push_back((i-0.5)/N);
    
    double vt {0.04}, dt {0.1};
    model_params::HestonParams heston_param {.S_0 = 100, .T = 10, .r = 0.0, .kappa = 1, .gamma = 0.3, .v_bar = vt, .v0 = 0.04, .rho = -0.9};
    double e = std::exp(-dt*heston_param.kappa);
    double m = e*(vt - heston_param.v_bar) + heston_param.v_bar;
    double var = heston_param.gamma*heston_param.gamma/heston_param.kappa*(heston_param.v_bar/2 * (1 - e)*(1-e) + vt*e*(1-e));
    double psi = var/(m*m);
    ASSERT_LT(psi, heston_mc::psi_c);

    for (auto u: u_vec)
        vt_vec.push_back(heston_mc::heston_mc_variance(vt, dt, heston_param,  u));
    double sample_mean =  std::accumulate(vt_vec.begin(), vt_vec.end(), 0.0) / vt_vec.size();

    double sum = 0.0;
    for (double x : vt_vec)
        sum += (x - sample_mean) * (x - sample_mean);
    double sampel_var = sum / vt_vec.size();
        
    
    EXPECT_NEAR(m, sample_mean, 1e-3*m);
    EXPECT_NEAR(var, sampel_var, 1e-2*var);
}

TEST(HestonMC, MomentMatchExponential){
    int N {100000};
    std::vector<double> u_vec, vt_vec{};
    for (int i=1; i<=N; ++i)
        u_vec.push_back((i-0.5)/N);
    
    double vt {0.0}, dt {0.1};
    model_params::HestonParams heston_param {.S_0 = 100, .T = 10, .r = 0.0, .kappa = 1, .gamma = 1, .v_bar = 0.04, .v0 = 0.04, .rho = -0.9};
    double e = std::exp(-dt*heston_param.kappa);
    double m = e*(vt - heston_param.v_bar) + heston_param.v_bar;
    double var = heston_param.gamma*heston_param.gamma/heston_param.kappa*(heston_param.v_bar/2 * (1 - e)*(1-e) + vt*e*(1-e));
    double psi = var/(m*m);
    ASSERT_GT(psi, heston_mc::psi_c);

    for (auto u: u_vec)
        vt_vec.push_back(heston_mc::heston_mc_variance(vt, dt, heston_param,  u));
    double sample_mean =  std::accumulate(vt_vec.begin(), vt_vec.end(), 0.0) / vt_vec.size();

    double sum = 0.0;
    for (double x : vt_vec)
        sum += (x - sample_mean) * (x - sample_mean);
    double sampel_var = sum / vt_vec.size();
        
    
    EXPECT_NEAR(m, sample_mean, 1e-3*m);
    EXPECT_NEAR(var, sampel_var, 1e-2*var);
}

TEST(HestonMC, NonNegative){
    int N {100000};
    std::vector<double> u_vec, vt_vec{};
    for (int i=1; i<=N; ++i)
        u_vec.push_back((i-0.5)/N);
    
    double vt {0.04}, dt {0.1};
    model_params::HestonParams heston_param {.S_0 = 100, .T = 10, .r = 0.0, .kappa = 1, .gamma = 1, .v_bar = 0.04, .v0 = 0.04, .rho = -0.9};
    std::vector<double> gamma_vec {0.3, 1.0};
    for (auto gamma:gamma_vec){
        heston_param.gamma = gamma;
        for (auto u: u_vec)
           EXPECT_GE( heston_mc::heston_mc_variance(vt, dt, heston_param,  u), 0.0) << "u = " << u;
    }
}

TEST(HestonMC, ExpectZeroBelowP){
    double vt {0.0}, dt {0.1};
    model_params::HestonParams heston_param {.S_0 = 100, .T = 10, .r = 0.0, .kappa = 1, .gamma = 1, .v_bar = 0.04, .v0 = 0.04, .rho = -0.9};
    double e = std::exp(-dt*heston_param.kappa);
    double m = e*(vt - heston_param.v_bar) + heston_param.v_bar;
    double var = heston_param.gamma*heston_param.gamma/heston_param.kappa*(heston_param.v_bar/2 * (1 - e)*(1-e) + vt*e*(1-e));
    double psi = var/(m*m);
    double p {(psi-1)/(psi+1)};

    int N = 10000;
    std::vector<double> u_vec;
    for (int i=1; i<=N*p; ++i)
        u_vec.push_back((i-0.5)/N);
    for (auto u: u_vec)
        EXPECT_EQ(heston_mc::heston_mc_variance(vt, dt, heston_param,  u), 0.0) << "u = " << u;
    EXPECT_GT(heston_mc::heston_mc_variance(vt, dt, heston_param, 0.99), 0.0);
}

TEST(HestonMC, BranchSelection){
    double dt {0.1};
    model_params::HestonParams heston_param {.S_0 = 100, .T = 10, .r = 0.0, .kappa = 1, .gamma = 1, .v_bar = 0.04, .v0 = 0.04, .rho = -0.9};

    // psi as a function of v_t for these params (Andersen moments)
    auto psi_of = [&](double vt){
        double e  = std::exp(-dt*heston_param.kappa);
        double g2 = heston_param.gamma*heston_param.gamma;
        double m  = heston_param.v_bar + (vt - heston_param.v_bar)*e;
        double s2 = vt*g2*e*(1-e)/heston_param.kappa
                  + heston_param.v_bar*g2*(1-e)*(1-e)/(2*heston_param.kappa);
        return s2/(m*m);
    };

    double u {0.1};

    // just below psi_c (~1.38): quadratic branch -> strictly positive at u = 0.1
    double vt_quad {0.07};
    ASSERT_LT(psi_of(vt_quad), heston_mc::psi_c);
    EXPECT_GT(heston_mc::heston_mc_variance(vt_quad, dt, heston_param, u), 0.0);

    // just above psi_c (~1.71, p ~ 0.26): exponential branch -> exactly 0 at u = 0.1 <= p
    double vt_exp {0.055};
    ASSERT_GT(psi_of(vt_exp), heston_mc::psi_c);
    EXPECT_EQ(heston_mc::heston_mc_variance(vt_exp, dt, heston_param, u), 0.0);
}
