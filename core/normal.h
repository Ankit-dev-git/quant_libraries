#pragma once
namespace normal {
	double randn();
	double norm_cdf(double x, double mu = 0, double sigma =1);
	double norm_pdf(double x, double mu = 0, double sigma =1);
}