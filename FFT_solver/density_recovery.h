#include <vector>
#include <functional>
#include <complex>
#include <FFT.h>

// Characteristic functions signature
using CF = std::function<Complex(double)>;

std::vector<double> RecoverDensity(CF cf, const std::vector<double>& x, int N = 64);