#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>
#include "FFT.h"
#include "density_recovery.h"
#include <vector>
#include <limits>

namespace{
    fourier::CF test_cf = [](double u)
    {
        Complex i(0.0, 1.0);
        double mu = 1, sigma = 0.2;
        return std::exp(i * mu * u - 0.5 * sigma * sigma * u * u);
    };
}

TEST(FFT, RecoverDensityMustFailN){
    std::vector<int> must_fail {5, 0, -8, -256, 12};
    for (int elem: must_fail)
        EXPECT_THROW(RecoverDensity(test_cf, 0.1, 0.0, elem), std::invalid_argument);
}

TEST(FFT, RecoverDensityMustPassN){
    std::vector<int> must_pass {1, 1024, 2};
    for (int elem: must_pass)
        EXPECT_NO_THROW(RecoverDensity(test_cf, 0.1, 0.0, elem));
}

TEST(FFT, RejectNonPowerN){
    Complex i (0.0, -1.0);
    std::vector<Complex> vec1(12, i), vec2(5, i);
    EXPECT_THROW(FFT(vec1), std::invalid_argument);
    EXPECT_THROW(FFT(vec2), std::invalid_argument);
}

TEST(FFT, RejectBadDu){
    std::vector<double> du {0.0, -0.7, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    for (double elem: du)
        EXPECT_THROW(RecoverDensity(test_cf, elem, 0.0, 64), std::invalid_argument);
}