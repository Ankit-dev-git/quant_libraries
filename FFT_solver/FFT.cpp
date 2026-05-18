#include "fft.h"
#include <cmath>

const double PI = 3.14159265358979323846;

void FFT(std::vector<Complex>& x)
{
    int N = x.size();

    if (N <= 1)
        return;

    std::vector<Complex> even(N / 2);
    std::vector<Complex> odd(N / 2);

    for (int i = 0; i < N / 2; ++i)
    {
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }

    FFT(even);
    FFT(odd);

    for (int k = 0; k < N / 2; ++k)
    {
        Complex twiddle =
            std::polar(1.0, -2.0 * PI * k / N);

        Complex t = twiddle * odd[k];

        x[k] = even[k] + t;
        x[k + N / 2] = even[k] - t;
    }
}