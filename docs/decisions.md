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