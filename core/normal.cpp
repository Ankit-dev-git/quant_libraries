#include "normal.h"
#include <random>
#include <cmath>

static thread_local std::mt19937 rng{ std::random_device{}() };

double normal::randn() {
    static thread_local std::normal_distribution<double> dist(0.0, 1.0);
    return dist(rng);
}

double normal::norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}