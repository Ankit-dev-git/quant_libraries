#pragma once
#include <vector>

enum class OptionType {call=1, put=-1};
namespace bs {
	double option_price(OptionType CP, double stock_price, double k, double sigma, double t, double T, double r);
	double option_delta(OptionType CP, double stock_price, double k, double sigma, double t, double T, double r);
	double d1(double S, double K, double sigma, double t, double T, double r);
}