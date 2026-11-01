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


## <17 Oct 2026> — usage of c4 in deciding the width of x-grid along with c2

- **Context:** the x-grid based on c2 was not able to cover all the non zero density points. Heston ReferencePrice at L=8/12 gave 5.8386 vs 5.785155 and the end-density test failed (fat tails).
- **Options:**
  1. Keep the grid based on c2 and use one large global L (L=16 worked) to ensure end point density is zero
  2. Incorporate c4 (kurtosis) so that grid adjusts automatically (Fang–Oosterlee: half-width L·√(c2 + √c4))
- **Decision:** Option 2. c4 < 0 is clamped to 0 in build_grid (estimator stays honest); non-finite c4 throws. c4 uses its own step h2 = 1e-2 (round-off ~16ε/h⁴ at the c2 step 1e-4 would swamp it).
- **Why:** Option 2 was selected because a fixed large L is a guess: it wastes resolution on thin-tailed models (BS) and may still be too small for the next fat-tailed model. c4 comes from each model's own CF, so the grid adapts automatically. Without adjusting, there was wrap around issue and density at higher x were non zero due to which call option price deviated.
- **Consequence:** Now the grid is adjusted by adjusting the effective sigma, which ensures that end point x have zero densities. BS is unchanged (c4 = 0); the Heston T=10 stress case converges (27.9090 vs 27.91). Cost: A coarser dx will require increasing N or decrease in U range.


## 17 Oct 2026 — Using put+parity route to price call options

- **Context:** After the c4 adjustment to the x-grid window, large call option values when the Feller condition is violated: huge density near S_T ~ 0, which wrapped around to large x-grid points.
- **Options:**
  1. Price the call directly
  2. Price the put, which integrates (K − eˣ)⁺, and get the call from C = P + S₀ − K·e^(−rT)
  3. Widen the grid or raise N until the direct call converges
- **Decision:** Option 2.
- **Why:** The call payoff (eˣ − K)⁺ grows exponentially at the right edge of the grid, exactly where the wrapped-around density error sits, so a small density error becomes a large price error. The put payoff is zero for x > ln K and bounded by K elsewhere, so the same error barely moves the put. Option 3 fights the eˣ growth with cost and doesn't remove it. The direct call is kept as a tail diagnostic.
- **Consequence:** Call pricing now depends on put–call parity, which in turn depends on the discounted stock price being a martingale under the CF. If a CF has a wrong drift or compensator, the call is wrong without warning, which is why martingale and reference tests matter. Calls are now as accurate as puts.


## 20 Oct 2026 — Dropping Feller enforcement while pricing through model characteristic function using FFT

- **Context:** The Feller condition is 2κv̄ ≥ γ². When it fails, the variance process can hit zero (the stock price never goes negative: S = eˣ). With Feller violated, at long T and large vol of vol, we see a large density mass near S_T ~ 0 when pricing through the CF with the FFT.
- **Options:**
  1. Throw when Feller is violated
  2. Warn but proceed
  3. Accept any valid parameters and enforce only the hard constraints: κ, v̄, γ > 0, v₀ ≥ 0, ρ ∈ (−1, 1).
- **Decision:** Option 3.
- **Why:** The CF is well defined even when Feller fails, since zero is a legitimate (reflecting) state for the variance. Calibrated equity parameters routinely violate Feller (steep skew needs a large γ), so enforcing it would block realistic fits. Feller violation mainly affects Monte Carlo, where a discretised variance step can go negative; that is handled in the MC scheme (truncation or Andersen QE), not by rejecting parameters.
- **Consequence:** FellerViolatedStillValid is a test for the Heston CF (and indirectly for Bates). The MC scheme must keep v ≥ 0 by construction.


## 20 Oct 2026 — Validate inputs only at the public boundary

- **Context:** Public functions take raw numbers; bad input (N not a power of 2, du ≤ 0, NaN params) would give unknown behaviour instead of failing loudly.
- **Options:**
  1. Validate in every function, internals included
  2. Validate at public entry points only; internals (fft_impl) in an anonymous namespace trust their inputs
  3. assert
- **Decision:** Option 2, throwing std::invalid_argument.
- **Why:** Option 1 re-checks the same inputs (e.g. N is a power of 2) at every recursion level of the FFT. Option 3 disappears in release builds, which CI runs. Internals in the anonymous namespace have internal linkage (callable only from that .cpp), so every call goes through the validating public wrapper and the checks can't be bypassed.
- **Consequence:** Every public function has a Bad*/Throws* test. Gap: params are validated in the CF builder, not at construction — planned fix: named constructor (static create() + private ctor) giving immutable, validated param types.


## 20 Oct 2026 — Returning cumulants along with characteristic function using CharacteristicFunctResult struct.

- **Context:** build_grid needs the centre and width of X_T, so a bare characteristic function would not have been enough and we would have had to calculate them in grid.
- **Options:**
  1. bare CF, with the grid builder computing cumulants numerically every time;
  2. a result struct bundling the CF with μ, σ, c₄ and the model name;
  3. the grid builder takes each model’s params directly.
- **Decision:** Option 2.
- **Why:** This ensured the grid builder only did the work it was meant to do, which is building a grid. Cumulants belong to the specific characteristic function and thus should be attached to it. Option 1 recomputes cumulants numerically on every call and throws away exact closed forms where they exist (BS, Merton). Option 3 couples build_grid to every model: each new model means editing the grid builder.
- **Consequence:** Every new model must fill all fields. A non-finite field makes build_grid throw, but a zero or wrong μ or c₄ is finite, so it silently gives a mis-centred or too-narrow grid; however this are guarded by TerminalMoments / AnalyticCumulants / AdditiveCumulants unit tests. The CF captures params by copy, so it goes stale if params change after building.

## 1 Nov 2026 — QE variance step as a pure function with the uniform passed in, not drawn inside

- **Context:** To compute the next Heston variance with Andersen QE, we need to design the step function. Used a pure function, which takes Heston params, vt, dt and one u value.
- **Options:**
  1. Both u and Z drawn inside the function (RNG inside).
  2. Function takes u for the exponential branch and Z for the quadratic branch as two separate inputs.
  3. Function takes one u from outside along with Heston params, dt and vt; then computes psi and, in the quadratic branch, Z = Φ⁻¹(u) using the inverse normal CDF (Acklam algorithm + one Halley step).
  4. Same as 3, but takes two uniforms and computes Z with the Box–Muller method.
- **Decision:** Option 3.
- **Why:** Option 1 cannot be tested against exact values from outside (only statistical tests). Options 2 and 3 are similar, but in option 2 the two branches use unrelated draws, so when psi crosses psi_c the output jumps between independent random numbers, which makes common-random-number Greeks noisy; option 3 drives both branches from the same u. Option 4 needs two uniforms per normal and the u → Z map is not monotone.
- **Consequence:** The function is pure: one u gives one v (reproducible, testable with a deterministic u-grid), and the u → Z map is monotone. Cost: Φ⁻¹ is an extra log/sqrt plus a Halley step (with norm_cdf) per call, so it is slower than drawing a normal directly; accuracy is ~machine precision after the Halley step.


## 1 Nov 2026 — BS, Heston, Merton and Bates parameter structs moved from FFT_solver/ to core/

- **Context:** All parameter structs lived in the recover_density library (FFT_solver/characteristic_func.h), so any library that wanted to use them had to link recover_density.
- **Options:**
  1. Create a new parameter struct inside the new library and use it directly there.
  2. Link recover_density into other libraries so they use the existing parameter struct for every model.
  3. Move all model parameter structs to the core library.
- **Decision:** Option 3.
- **Why:** Option 1 gives flexibility and no dependency on the Fourier library, but leads to duplicate structs for the same models across the pricing library. Option 2 would link recover_density into libraries that don't use any of its functionality. Option 3 gives every model one struct in core. Core holds shared basics (normal distribution functions, BS analytics) and is linked by almost every library anyway. Aliases in characteristic_func.h keep fourier::HestonParams etc. valid, so existing code and tests are unchanged.
- **Consequence:** One struct per model keeps the pricing library clean; recover_density now links core PUBLIC because its public header includes a core header. Cost: a shared struct must hold every field any user needs (e.g. S_0 and T are unused by the MC variance step), so it can grow; the structs are still plain aggregates validated in the CF builders 


## 1 Nov 2026 — QE vs exact non-central χ² sampling

- **Context:** We can sample v_{t+Δ} exactly from the non-central chi-square distribution (Grzelak's course), but we chose the Andersen QE approximation instead.
- **Options:**
  1. Exact non-central chi-square sampling of v_{t+Δ} (Poisson mixture of chi-squares).
  2. Andersen QE: match the exact conditional mean and variance with a simple distribution.
  3. Euler with truncation/reflection at zero.
- **Decision:** Option 2.
- **Why:** QE is much faster than exact sampling (a few closed-form operations on one uniform vs a Poisson draw plus a gamma/chi-square draw with varying shape per step) and is the industry standard. QE maps one uniform to v through an inverse CDF, so a small parameter bump moves v smoothly under common random numbers; the exact sampler's Poisson count jumps discretely when a parameter is bumped, which makes bump-and-reprice Greeks noisy. Euler is biased over a finite step and can go negative, so it needs ad-hoc fixes.
- **Consequence:** Only the first two moments are exact, so there is a small discretisation bias.