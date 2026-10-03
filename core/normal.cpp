#include "normal.h"
#include <random>
#include <cmath>
#include<stdexcept>

static thread_local std::mt19937 rng{ std::random_device{}() };

double normal::randn() {
    static thread_local std::normal_distribution<double> dist(0.0, 1.0);
    return dist(rng);
}

double normal::norm_cdf(double x, double mu, double sigma) {
    if (sigma <= 0)
        throw std::invalid_argument("norm cdf: sigma must be positive");
    double z {(x-mu)/sigma};
    return 0.5 * std::erfc(-z / std::sqrt(2.0));
}

double normal::norm_pdf(double x, double mu, double sigma){
    if (sigma <= 0)
        throw std::invalid_argument("norm pdf: sigma must be positive");
    double z {(x-mu)/sigma};
    return 1/(std::sqrt(2.0*std::numbers::pi)*sigma)*std::exp(-z*z/2);
}