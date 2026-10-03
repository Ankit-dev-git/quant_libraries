# quant_libraries

![CI](https://github.com/Ankit-dev-git/quant_libraries/actions/workflows/ci.yml/badge.svg)

A C++20 derivatives pricing library built from scratch: Monte Carlo simulation, discrete delta hedging, and Fourier methods (characteristic-function inversion via FFT). Built with CMake, tested with GoogleTest, and continuously integrated on Linux via GitHub Actions.

---

## Projects

### 1. Dynamic Delta Hedging (`Dynamic_hedging_BS/`)

Simulates the discrete delta hedging of a European option under Black–Scholes and analyses the resulting hedging P&L.

- **Model:** Geometric Brownian Motion paths; Black–Scholes price and delta recomputed at each rebalancing date
- **What it measures:** distribution of hedging P&L across paths
- **Performance:** delta-hedging simulation parallelised with a custom `parallel_for` thread pool
  - 50,000 paths × 1,000 steps, Release build (GCC 15, WSL2): **7.5 s single-threaded → 1.7 s with 8 threads on 4 physical cores (4.3×)**
  - Median of 5 runs; timings cover the per-path hedging loop only (GBM path generation runs beforehand and is not timed). Single- and multi-threaded P&L outputs are identical. The gain above 4× is most likely from hyperthreading (run order and thermal effects were checked and ruled out).

### 2. Fourier density recovery (`FFT_solver/`)

Recovers a probability density from its characteristic function using a radix-2 Cooley–Tukey FFT with a trapezoidal boundary correction, checked against the exact normal density.
*In progress:* CMake integration and GoogleTest checks. Next: COS option pricing with Black–Scholes and Heston characteristic functions.

---

## Build and test

Requirements: a C++20 compiler (tested with GCC 13 in CI on Ubuntu 24.04, and GCC 15), CMake ≥ 3.25, Ninja.

```bash
git clone https://github.com/Ankit-dev-git/quant_libraries.git
cd quant_libraries

cmake -S . -B build-Release -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-Release
ctest --test-dir build-Release --output-on-failure
```

For a debug build, replace `Release` with `Debug` in the three commands.

Pass `-DBUILD_TESTING=OFF` at configure time to skip the tests (GoogleTest is then not downloaded). Otherwise GoogleTest is fetched automatically at configure time; no manual installation needed.

## Run

```bash
./build-Release/Dynamic_hedging_BS/dynamic_hedging
```

Prints single- and multi-threaded timings. Output CSV files are written to the current working directory.

---

## Repository structure

```
quant_libraries/
├── CMakeLists.txt                  # top-level build: settings, warnings, GoogleTest
├── .github/workflows/ci.yml        # build + test on every push (warnings are errors)
├── Dynamic_hedging_BS/
│   ├── CMakeLists.txt              # 'hedging' library + 'dynamic_hedging' executable
│   ├── Dynamic_hedging_BS.cpp      # main(): runs and times the simulation
│   ├── BS_price.{h,cpp}            # Black–Scholes price, delta, d1
│   ├── generate_paths.{h,cpp}      # GBM path simulation
│   ├── delta_hedge.{h,cpp}         # discrete hedging and P&L
│   ├── parallel_for.{h,cpp}        # thread pool / parallel loop
│   ├── GridCreation.{h,cpp}        # flattened 2-D grid for path storage
│   ├── DataContainerCreation.{h,cpp}
│   ├── normal.{h,cpp}              # normal CDF / PDF
│   └── tests/                      # GoogleTest unit tests
└── FFT_solver/                     # not yet in the build
    ├── FFT.{h,cpp}                 # recursive radix-2 FFT
    └── FFT_solver.cpp              # CF → density recovery + main()
```

---

## Testing

Unit tests use GoogleTest and run in CI on every push. Current coverage:

- **Known value:** Black–Scholes call with S = K = 100, r = 5%, σ = 20%, T = 1 → 10.4506 (tolerance 1e-4)
- **Put–call parity:** C − P = S − K·e^(−rT) to 1e-10
- **At expiry:** price equals payoff and delta is the step function, for calls and puts in, at and out of the money

CI builds in Release with compiler warnings treated as errors (`-Wall -Wextra -Wpedantic`).

Tests are added alongside every new model.