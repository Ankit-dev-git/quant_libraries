#include "normal.h"
#include <random>
#include <cmath>
#include <numbers>
#include<stdexcept>

static thread_local std::mt19937 rng{ std::random_device{}() };

namespace {
    // Acklam (2003) rational-approximation coefficients for the inverse normal CDF
    constexpr double a[] = {-3.969683028665376e+01,  2.209460984245205e+02, -2.759285104469687e+02,
                             1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00};
    constexpr double b[] = {-5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,
                             6.680131188771972e+01, -1.328068155288572e+01};
    constexpr double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                            -2.549732539343734e+00,  4.374664141464968e+00,  2.938163982698783e+00};
    constexpr double d[] = { 7.784695709041462e-03,  3.224671290700398e-01,  2.445134137142996e+00,
                             3.754408661907416e+00};
    constexpr double p_low = 0.02425;

    // tail rational function in q = sqrt(-2 ln p); used for both tails (upper = mirrored)
    double inv_cdf_tail(double q){
        return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5])
             / ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
    }
}

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

double normal::norm_inv_cdf(double p){
    if (!(p > 0.0 && p < 1.0))      // also rejects NaN
        throw std::invalid_argument("norm_inv_cdf: p must be in the open interval (0,1)");

    double x;
    if (p < p_low){                                   // lower tail
        x = inv_cdf_tail(std::sqrt(-2.0 * std::log(p)));
    }
    else if (p <= 1.0 - p_low){                       // central region
        double q = p - 0.5;
        double r = q * q;
        x = (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q
          / (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0);
    }
    else{                                             // upper tail, mirror of the lower
        x = -inv_cdf_tail(std::sqrt(-2.0 * std::log(1.0 - p)));
    }

    // one Halley refinement step: ~1e-9 -> ~machine precision
    double e = norm_cdf(x) - p;
    double u = e * std::sqrt(2.0 * std::numbers::pi) * std::exp(0.5 * x * x);
    x = x - u / (1.0 + 0.5 * x * u);
    return x;
}