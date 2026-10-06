#pragma once

#include <vector>
#include <complex>
#include <numbers>

using Complex = std::complex<double>;
namespace fourier{
    constexpr double PI = std::numbers::pi;
    using CF = std::function<Complex(double)>;
}

void FFT(std::vector<Complex>& x);