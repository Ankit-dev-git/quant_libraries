# FFT_solver (namespace fourier)

## 4 Oct 2026 — dx computed inside RecoverDensity, not passed in

- **Context:** dx and du are tied by dx·du = 2π/N. If a caller passes both, they can pass an inconsistent pair.
- **Options:**
  1. Build dx (and the x-grid) in `main` and pass it to `RecoverDensity`.
  2. Pass dx in, but validate it (throw if it doesn't equal 2π/(N·du)).
  3. Compute dx inside `RecoverDensity` from N and du, and return the grid in `DensityResult`.
- **Decision:** Option 3.
- **Why:** Once N and du are fixed, dx has exactly one valid value, so a dx parameter can only be right or silently wrong. Tested on the normal density: correct dx → max error 1.4e-16; dx 1% off → 2.6e-2, with no warning. Option 2 catches that, but a parameter with only one allowed value is redundant.
- **Consequence:** Callers use the returned grid, so the pricer and `main` can't drift from it. Cost: the caller can't choose dx directly; resolution is controlled through du and N (dx = 2π/(N·du)).


## 6 Oct 2026 — Model parameters passed as a struct

- **Context:** How to pass parameters to the characteristic function `charac_func_BS`.
- **Options:**
  1. Pass each parameter directly to the characteristic function.
  2. Pass all parameters as a struct.
- **Decision:** Option 2.
- **Why:** Passing parameters individually risks passing them in the wrong order. It would also mean creating several separate variables in `main`, since the parameters are used downstream too.
- **Consequence:** No mismatch from argument order, and the struct keeps all parameters in one place. Cost: each model needs its own struct type. 


## 6 Oct 2026 — CF factory over inheritance

- **Context:** Whether `RecoverDensity` takes a `fourier::CF` or a virtual base class as its characteristic-function argument.
- **Options:**
  1. Create a virtual base class for the characteristic function, with each model's characteristic function derived from it.
  2. Create a factory: each model's characteristic function is built by its own function and passed to `RecoverDensity` as a `fourier::CF`.
- **Decision:** Option 2.
- **Why:** A virtual base class would need a derived class for each model, with its parameters as members. With the factory, each model needs only one function, and `RecoverDensity` can take any callable as `fourier::CF`.
- **Consequence:** Each model's characteristic function has full control over how it is built. Cost: there are several functions to know and call, since each model's characteristic function is built by a different function. Also, `mu` and `sigma_XT` are captured by copy, so once the CF is created, changing any parameter means recreating the CF.


## 6 Oct 2026 — dx read from the grid in the pricer

- **Context:** The pricer can recreate the x-grid from du, N, x_min and u_span, or it can take the grid from `RecoverDensity`, which already builds it to compute the density at those points with the FFT.
- **Options:**
  1. Use du, N, u_span and x_min (and so dx) directly, and recreate the grid in the pricer.
  2. Return the grid in the `RecoverDensity` result. The pricer reads dx directly from the grid in the `DensityResult` struct returned by the density function.
- **Decision:** Option 2.
- **Why:** Recreating the x-grid could mean pricing on a different grid from the one the density was computed on. This can happen if the density function changes the grid internally (e.g. shifts x_min) and returns a modified grid. A dx 1% off between pricer and density gave an error of 2.6e-2, with no warning.
- **Consequence:** The density function has full control over the x-grid used in the pricer. Cost: the pricer takes x[1] − x[0] as dx and assumes dx is constant. Also, if the grid has 0 or 1 points, the pricer cannot compute dx and throws an error.


## 6 Oct 2026 — Model-independent pricer

- **Context:** The pricer can take each model's characteristic function and price from it, or it can take the terminal density and parameters such as K, r and T.
- **Options:**
  1. A model-dependent pricer for each model, taking that model's characteristic function as input.
  2. One pricer that takes the terminal density and grid from `DensityResult` (model-dependent) plus market inputs (K, r, T).
- **Decision:** Option 2.
- **Why:** The pricer has been tested against the BS closed form, so any error in another model can be attributed to its density or characteristic function, not the pricer. One pricer also keeps the architecture clean, and the tests built on it don't need to change for every model.
- **Consequence:** Adding or changing a model does not require a new pricer. Cost: because the pricer only sees the density of S_T at maturity, path-dependent derivatives cannot be priced with it.


## <date> — <decision in a few words>

- **Context:** <the problem or constraint that forced a choice — not the answer>
- **Options:**
  1. <alternative A>
  2. <alternative B>
- **Decision:** Option <n>.
- **Why:** <the one main argument>. <evidence: a measured number, test result or concrete failure case>.
- **Consequence:** <what this makes easier or safer>. Cost: <what it makes harder or rules out>.