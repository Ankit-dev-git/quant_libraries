#include "BS_price.h"
#include "normal.h"
#include <cmath>

double bs::option_price(OptionType CP, double stock_price, double k, double sigma, double t, double T, double r)
{	
	double D1{ 0 }, value{ 0 };
	D1 = d1(stock_price, k, sigma, t, T, r);
	double d2{ D1 - sigma * std::sqrt(T - t) };
	if (CP == OptionType::call)
		value = stock_price * normal::norm_cdf(D1) - k * normal::norm_cdf(d2) * std::exp(-r * (T - t));
	else if (CP == OptionType::put)
		value = k * normal::norm_cdf(-d2) * std::exp(-r * (T - t)) - stock_price * normal::norm_cdf(-D1);

	return value;
}

double bs::option_delta(OptionType CP, double stock_price, double k, double sigma, double t, double T, double r)
{
	double D1{ 0 }, delta{ 0 } ;
	D1 = d1(stock_price, k, sigma, t, T, r);
	if (CP == OptionType::call)
		delta = normal::norm_cdf(D1);
	else if (CP == OptionType::put)
		delta = normal::norm_cdf(D1) - 1.0;
	return delta;
}

double bs::d1(double stock_price, double k, double sigma, double t, double T, double r)
{
	return (std::log(stock_price / k) + (r + sigma * sigma / 2) * (T - t)) / (sigma * std::sqrt(T - t));
}
