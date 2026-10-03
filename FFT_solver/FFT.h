#pragma once

#include <vector>
#include <complex>
#include <numbers>

using Complex = std::complex<double>;
constexpr double PI = std::numbers::pi;

void FFT(std::vector<Complex>& x);