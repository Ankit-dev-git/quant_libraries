#include <gtest/gtest.h>
#include "grid_builder.h"
#include "characteristic_func.h"
#include "density_recovery.h"
#include <limits>
#include <vector>
#include <cmath>

TEST(GridBuilder, LnKOnGrid){
    fourier::BlackScholes_params bs_param {.S_0 = 100.0, .r = 0.05, .sigma = 0.2, .T = 1};
    fourier::CharacteristicFuncResult BS_charac = fourier::charac_func_BS(bs_param);
    int N = 64;
    double K = 100.0, L = 4;
    fourier::GridParams grid_param =  fourier::build_grid(BS_charac, K, L, N);
    double dx = (2*fourier::PI)/(N*grid_param.du);
    fourier::DensityResult density = fourier::RecoverDensity(BS_charac.charac_func, grid_param.du, grid_param.x_min, N);

    double K_ITM = 20;
    fourier::GridParams grid_param_ITM =  fourier::build_grid(BS_charac, K_ITM, L, N);

    double K_OTM = 400;
    fourier::GridParams grid_param_OTM =  fourier::build_grid(BS_charac, K_OTM, L, N);

    ASSERT_GE(grid_param.n_dx_xmin,0);
    ASSERT_LE(grid_param.n_dx_xmin, N-1);
    EXPECT_NEAR(grid_param.x_min + grid_param.n_dx_xmin*dx, std::log(K), 1e-12);
    EXPECT_NEAR(density.x[grid_param.n_dx_xmin], std::log(K), 1e-12);

    EXPECT_LT(grid_param_ITM.n_dx_xmin, 0);
    EXPECT_LE(std::log(K_ITM), grid_param_ITM.x_min);

    EXPECT_GT(grid_param_OTM.n_dx_xmin, N-1);
    EXPECT_GE(std::log(K_OTM), grid_param_OTM.x_min + (N-1)*dx);

}

    
TEST(GridBuilder, CenteringGridOnMu){
    fourier::BlackScholes_params bs_param {.S_0 = 100.0, .r = 0.05, .sigma = 0.2, .T = 1};
    fourier::CharacteristicFuncResult BS_charac = fourier::charac_func_BS(bs_param);
    int N = 64;
    double K = 100.0, L = 4;
    fourier::GridParams grid_param1 =  fourier::build_grid(BS_charac, K, L, N);
    double dx = (2*fourier::PI)/(N*grid_param1.du);
    fourier::DensityResult density1 = fourier::RecoverDensity(BS_charac.charac_func, grid_param1.du, grid_param1.x_min, N);

    double K_ITM = 20;
    fourier::GridParams grid_param_ITM =  fourier::build_grid(BS_charac, K_ITM, L, N);
    fourier::DensityResult density_ITM = fourier::RecoverDensity(BS_charac.charac_func, grid_param_ITM.du, grid_param_ITM.x_min, N);

    double K_OTM = 400;
    fourier::GridParams grid_param_OTM =  fourier::build_grid(BS_charac, K_OTM, L, N);
    fourier::DensityResult density_OTM = fourier::RecoverDensity(BS_charac.charac_func, grid_param_OTM.du, grid_param_OTM.x_min, N);

    EXPECT_LE(std::abs(grid_param1.x_min + N/2*dx - BS_charac.mu_terminal), 0.5*dx + 1e-12);
    EXPECT_LE(std::abs(density1.x[N/2] - BS_charac.mu_terminal), 0.5*dx + 1e-12);

    EXPECT_LE(std::abs(grid_param_ITM.x_min + N/2*dx - BS_charac.mu_terminal), 0.5*dx + 1e-12);
    EXPECT_LE(std::abs(density_ITM.x[N/2] - BS_charac.mu_terminal), 0.5*dx + 1e-12);

    EXPECT_LE(std::abs(grid_param_OTM.x_min + N/2*dx - BS_charac.mu_terminal), 0.5*dx + 1e-12);
    EXPECT_LE(std::abs(density_OTM.x[N/2] - BS_charac.mu_terminal), 0.5*dx + 1e-12);
}        

TEST(GridBuilder, BadN){
    fourier::BlackScholes_params bs_param {.S_0 = 100.0, .r = 0.05, .sigma = 0.2, .T = 1};
    fourier::CharacteristicFuncResult BS_charac = fourier::charac_func_BS(bs_param);
    std::vector<int> vec{0, -5};
    double K = 100.0, L = 4;
    for (auto& elem: vec)
        EXPECT_THROW(fourier::build_grid(BS_charac, K, L, elem), std::invalid_argument) << "N = " << elem;
}

TEST(GridBuilder, BadK){
    fourier::BlackScholes_params bs_param {.S_0 = 100.0, .r = 0.05, .sigma = 0.2, .T = 1};
    fourier::CharacteristicFuncResult BS_charac = fourier::charac_func_BS(bs_param);
    std::vector<double> vec{0, -5, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    double L = 4;
    for (auto& elem: vec)
        EXPECT_THROW(fourier::build_grid(BS_charac, elem, L), std::invalid_argument) << "K = " << elem;
}


TEST(GridBuilder, BadL){
    fourier::BlackScholes_params bs_param {.S_0 = 100.0, .r = 0.05, .sigma = 0.2, .T = 1};
    fourier::CharacteristicFuncResult BS_charac = fourier::charac_func_BS(bs_param);
    std::vector<double> vec{0, -5, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    double K = 100;
    for (auto& elem: vec)
        EXPECT_THROW(fourier::build_grid(BS_charac, K, elem), std::invalid_argument) << "L = " << elem;
}

TEST(GridBuilder, BadMu){
    fourier::BlackScholes_params bs_param {.S_0 = 100.0, .r = 0.05, .sigma = 0.2, .T = 1};
    fourier::CharacteristicFuncResult BS_charac = fourier::charac_func_BS(bs_param);
    std::vector<double> vec{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    double K = 100, L = 4;
    for (auto& elem: vec){
        BS_charac.mu_terminal = elem;
        EXPECT_THROW(fourier::build_grid(BS_charac, K, L), std::invalid_argument) << "Mu_terminal = " << elem;
    }
}

TEST(GridBuilder, BadSigma){
    fourier::BlackScholes_params bs_param {.S_0 = 100.0, .r = 0.05, .sigma = 0.2, .T = 1};
    fourier::CharacteristicFuncResult BS_charac = fourier::charac_func_BS(bs_param);
    std::vector<double> vec{0, -5, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    double K = 100, L = 4;
    for (auto& elem: vec){
        BS_charac.sigma_terminal = elem;
        EXPECT_THROW(fourier::build_grid(BS_charac, K, L), std::invalid_argument) << "Sigma_terminal = " << elem;
    }
}