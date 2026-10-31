# quant_libraries

![CI](https://github.com/Ankit-dev-git/quant_libraries/actions/workflows/ci.yml/badge.svg)

A C++20 derivatives pricing library built from scratch: Monte Carlo simulation, discrete delta hedging, and Fourier pricing (characteristic function → density via FFT) under Black–Scholes, Heston, Merton and Bates. Built with CMake, tested with GoogleTest, and continuously integrated on Linux via GitHub Actions.

**Headline results**
- Fourier pricer: O(Δx²) convergence with log-strike grid alignment, Richardson → O(Δx⁴) (**1.7e-8 error vs Black–Scholes at N = 1024**); Heston T = 10 stress case **27.9090 vs reference 27.91**; Merton **4.8e-4 vs the Merton (1976) series**.
- Monte Carlo delta hedging: **4.3× speed-up** (8 threads, 4 physical cores) with identical output.

---

## Projects

### 1. Dynamic Delta Hedging (`Dynamic_hedging_BS/`)

Simulates the discrete delta hedging of a European option under Black–Scholes and analyses the resulting hedging P&L.

- **Model:** Geometric Brownian Motion paths; Black–Scholes price and delta recomputed at each rebalancing date
- **What it measures:** distribution of hedging P&L across paths
- **Performance:** delta-hedging simulation parallelised with a custom `parallel_for` thread pool
  - 50,000 paths × 1,000 steps, Release build (GCC 15, WSL2): **7.5 s single-threaded → 1.7 s with 8 threads on 4 physical cores (4.3×)**
  - Median of 5 runs; timings cover the per-path hedging loop only (GBM path generation runs beforehand and is not timed). Single- and multi-threaded P&L outputs are identical. The gain above 4× is most likely from hyperthreading (run order and thermal effects were checked and ruled out).

### 2. Fourier pricer (`FFT_solver/`)

Prices European options by recovering the risk-neutral density of log S_T from the model's characteristic function with a radix-2 Cooley–Tukey FFT (trapezoidal boundary correction), then integrating the payoff against it.

- **Models:** Black–Scholes, Heston ("little trap" CF, stable for long maturities), Merton jump-diffusion, Bates (Heston × Merton jump factor)
- **Grid:** truncation window from the Fang–Oosterlee cumulant rule L·√(c₂ + √c₄), with c₂, c₄ computed numerically from the CF; centred on the mean, log-strike aligned to a grid node
- **Calls via puts + put–call parity**, because the eˣ factor amplifies density error in the right tail
- **Accuracy:**
  - Black–Scholes: O(Δx²); Richardson extrapolation → O(Δx⁴), **1.7e-8 at N = 1024**
  - Heston, T = 10, ρ = −0.9: **27.9013 / 27.9075 / 27.9090 at N = 1024 / 2048 / 4096** vs reference 27.91; step ratio 4.00; Richardson 27.9095
  - Merton: **4.8e-4 vs the Merton (1976) series at N = 1024**
  - Bates: reduces to Heston at λ = 0 and to Merton at small vol-of-vol (tested)
- **Input validation** at the public boundary (power-of-2 N, grid spacing, model parameters); internals in anonymous namespaces

*Next:* Heston Monte Carlo (Andersen QE) cross-checked against the Fourier price, then Heston calibration.

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
├── docs/decisions.md               # design decisions: context / options / decision / why / consequence
├── core/
│   ├── CMakeLists.txt
│   ├── BS_price.{h,cpp}            # Black–Scholes price, delta, d1
│   ├── normal.{h,cpp}              # normal CDF / PDF
│   └── tests/                      # GoogleTest unit tests
├── Dynamic_hedging_BS/
│   ├── CMakeLists.txt              # 'hedging' library + 'dynamic_hedging' executable
│   ├── Dynamic_hedging_BS.cpp      # main(): runs and times the simulation
│   ├── generate_paths.{h,cpp}      # GBM path simulation
│   ├── delta_hedge.{h,cpp}         # discrete hedging and P&L
│   ├── parallel_for.{h,cpp}        # thread pool / parallel loop
│   ├── GridCreation.{h,cpp}        # flattened 2-D grid for path storage
│   ├── DataContainerCreation.{h,cpp}
│   └── tests/                      # GoogleTest unit tests
└── FFT_solver/
    ├── CMakeLists.txt
    ├── FFT.{h,cpp}                 # recursive radix-2 FFT (validating public wrapper)
    ├── characteristic_func.{h,cpp} # CFs: Black–Scholes, Heston, Merton, Bates; numerical cumulants
    ├── grid_builder.{h,cpp}        # cumulant-based truncation window, log-strike-aligned grid
    ├── density_recovery.{h,cpp}    # CF → density via FFT
    ├── fourier_pricer.{h,cpp}      # payoff integration against the density (put + parity for calls)
    ├── run_density_recovery.cpp    # demo executable
    └── tests/                      # CF, grid and size tests per model (BS, Heston, Merton, Bates)
```

---

## Testing

Unit tests use GoogleTest and run in CI on every push. Current coverage:

**Black–Scholes analytics**
- **Known value:** Black–Scholes call with S = K = 100, r = 5%, σ = 20%, T = 1 → 10.4506 (tolerance 1e-4)
- **Put–call parity:** C − P = S − K·e^(−rT) to 1e-10
- **At expiry:** price equals payoff and delta is the step function, for calls and puts in, at and out of the money

**Fourier pricer**
- **CF sanity (parameterised over models):** φ(0) = 1, |φ| ≤ 1, conjugate symmetry, decay
- **Reference prices:** Merton (1976) series, Heston long-dated stress case
- **Limits:** Heston → Black–Scholes; Merton at λ = 0 → Black–Scholes; Bates at λ = 0 → Heston; Bates → Merton at small vol-of-vol
- **Cumulants:** numerical c₂, c₄ vs analytic (Merton, Bates additive) and a synthetic polynomial CF
- **Grid & FFT:** FFT rejects non-power-of-2 sizes; log-strike lands on a grid node; truncation window widens with c₄, clamps c₄ < 0, throws on non-finite cumulants
- **Robustness:** Feller-violating Heston parameters price without throwing; invalid inputs throw

CI builds in Release with compiler warnings treated as errors (`-Wall -Wextra -Wpedantic -Wnon-virtual-dtor`).

Tests are added alongside every new model.
