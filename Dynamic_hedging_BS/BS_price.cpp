#include "BS_price.h"
#include "normal.h"
#include <cmath>
#include <algorithm>

double bs::option_price(OptionType CP, double stock_price, double k, double sigma, double t, double T, double r)
{	
	int unused = 0;
	double value{ 0 };
	if(T-t < 1e-12){
		if (CP == OptionType::call)
			value = std::max(stock_price - k,0.0);
		else if (CP == OptionType::put)
			value = std::max(k - stock_price, 0.0);
	}	
	else
	{
		double D1{ 0 };
		D1 = d1(stock_price, k, sigma, t, T, r);
		double d2{ D1 - sigma * std::sqrt(T - t) };
		if (CP == OptionType::call)
			value = stock_price * normal::norm_cdf(D1) - k * normal::norm_cdf(d2) * std::exp(-r * (T - t));
		else if (CP == OptionType::put)
			value = k * normal::norm_cdf(-d2) * std::exp(-r * (T - t)) - stock_price * normal::norm_cdf(-D1);
	}

	return value;
}

double bs::option_delta(OptionType CP, double stock_price, double k, double sigma, double t, double T, double r)
{
	double delta{ 0 };
	if(T-t < 1e-12){
		if (CP == OptionType::call)
			delta = stock_price>k ? 1.0:0.0;
		else if (CP == OptionType::put)
			delta = stock_price<k ? -1.0:0.0;
	}	
	else{
		double D1{ 0 };
		D1 = d1(stock_price, k, sigma, t, T, r);
		if (CP == OptionType::call)
			delta = normal::norm_cdf(D1);
		else if (CP == OptionType::put)
			delta = normal::norm_cdf(D1) - 1.0;
	}
	return delta;
}

double bs::d1(double stock_price, double k, double sigma, double t, double T, double r)
{
	return (std::log(stock_price / k) + (r + sigma * sigma / 2) * (T - t)) / (sigma * std::sqrt(T - t));
}
