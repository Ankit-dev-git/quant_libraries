# pricing_lib

![CI](https://github.com/<your-username>/pricing_lib/actions/workflows/ci.yml/badge.svg)

A C++20 derivatives pricing library built from scratch: Monte Carlo simulation, discrete delta hedging, and Fourier-based option pricing. Built with CMake, tested with GoogleTest, and continuously integrated on Linux via GitHub Actions.

---

## Projects

### 1. Dynamic Delta Hedging (`Dynamic_hedging_BS/`)

Simulates the discrete delta hedging of a European option under Black–Scholes and analyses the resulting hedging P&L.

- **Model:** Geometric Brownian Motion paths; Black–Scholes price and delta recomputed at each rebalancing date
- **What it measures:** distribution of hedging P&L across thousands of paths and time steps, and how hedging error depends on rebalancing frequency
- **Performance:** multi-threaded path simulation using a custom `parallel_for`
  - 50,000 paths: **9,181 ms single-threaded → 2,563 ms on 8 threads**
  - **3.6× speedup on 4 physical cores (~90% parallel efficiency)**

### 2. FFT Option Pricer (`FFT_solver/`)

Carr–Madan FFT pricing of European options using Simpson's rule quadrature, validated against the Black–Scholes closed form.

---

## Build and test

Requirements: a C++20 compiler (GCC ≥ 11 or Clang ≥ 14), CMake ≥ 3.20, Ninja (optional).

```bash
git clone https://github.com/<your-username>/pricing_lib.git
cd pricing_lib
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

GoogleTest is fetched automatically at configure time — no manual installation needed.

---

## Repository structure

```
pricing_lib/
├── CMakeLists.txt              # top-level build
├── .github/workflows/ci.yml    # build + test on every push
├── Dynamic_hedging_BS/
│   ├── CMakeLists.txt          # 'hedging' library + 'dynamic_hedging' executable
│   ├── BS_price.{h,cpp}        # Black–Scholes price, delta, d1
│   ├── generate_paths.{h,cpp}  # GBM path simulation
│   ├── delta_hedge.{h,cpp}     # discrete hedging and P&L
│   ├── parallel_for.{h,cpp}    # multi-threaded loop execution
│   ├── normal.{h,cpp}          # normal CDF / PDF
│   └── tests/                  # GoogleTest unit tests
└── FFT_solver/
    └── FFT.{h,cpp}             # Carr–Madan FFT pricer
```

---

## Testing

Unit tests use GoogleTest and run in CI on every push. Current coverage:

- Black–Scholes call price against a known reference value
- Put–call parity

Tests are added alongside every new model.

---

## Roadmap

- [ ] Heston, Merton jump-diffusion and Bates models via a shared characteristic-function interface (CRTP), priced with both COS and Carr–Madan
- [ ] Heston Monte Carlo, cross-checked against Fourier prices
- [ ] Barrier and Asian options in the Monte Carlo engine
- [ ] Yield curve bootstrapping (OIS)
- [ ] Algorithmic adjoint differentiation (AAD) for Greeks
