#include "FFT.h"
#include <cmath>
#include <vector>
#include <stdexcept>

// Creating anonymous namespace as this function is called only in this cpp. This is similar to static double fft_imp(...){...};
namespace
{
    void fft_impl(std::vector<Complex>& x){
        int N = x.size();
        if (N<=1)
            return;

        std::vector<Complex> even(N / 2);
        std::vector<Complex> odd(N / 2);

        for (int i = 0; i < N / 2; ++i)
        {
            even[i] = x[2 * i];
            odd[i] = x[2 * i + 1];
        }

        fft_impl(even);
        fft_impl(odd);

        for (int k = 0; k < N / 2; ++k)
        {
            Complex twiddle =
                std::polar(1.0, -2.0 * fourier::PI * k / N);

            Complex t = twiddle * odd[k];

            x[k] = even[k] + t;
            x[k + N / 2] = even[k] - t;
        }
    }
} 


void FFT(std::vector<Complex>& x)
{
    int N = x.size();
// Adding the check outside recursion so that it is not checked under every recursive call. Using private wrapper function fft_imp to achieve the same    
    if (N>0 && std::floor(std::log2(N))==std::log2(N)){
        if (N==1)
            return;
        fft_impl(x);
    }
    else
        throw std::invalid_argument("N should be a power of 2 where power is whole number");
}