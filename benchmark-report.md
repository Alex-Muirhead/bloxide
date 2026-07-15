# bloxide performance history: upstream/master -> HEAD

Per-commit benchmark of `bloxide examples/anderson.yaml` (`--release`) for every commit
on `master` since `upstream/master`, plus the current working tree (scan-iterator
reverted to a for-loop).

**Method:** each commit built in an isolated git worktree with a shared `target/`
(rustc 1.89.0), binaries stashed and benchmarked in a single interleaved
`hyperfine -w 3 -m 100` run. Machine: AMD Ryzen Threadripper PRO 3955WX, Linux
(Pop!_OS, kernel 7.0.11), hyperfine 1.20.0. All builds produce identical physical
results to printed precision (skin friction 347.52980 N/m^2, heat transfer
61.86821 W/cm^2, adiabatic wall temp 5035.86311 K); only the Newton iteration count
varies slightly (16-18) from floating-point reassociation.

## Results

| #  | Commit    | Subject                              | Mean +/- sigma    | vs upstream | vs previous |
|----|-----------|--------------------------------------|-------------------|-------------|-------------|
| 01 | `214fb25` | upstream/master (baseline)           | 37.43 +/- 0.41 ms | 1.00x       |  -          |
| 06 | `d9f262e` | Migrate over to num-dual             | 16.74 +/- 0.20 ms | 2.24x       | 2.23x       |
| 08 | `f0bf71e` | Use dual numbers for derivatives     | 17.76 +/- 0.35 ms | 2.11x       | 0.94x       |
| 10 | `0286b68` | Try moving to doing 2 derivs at once | 13.73 +/- 0.33 ms | 2.73x       | 1.30x       |
| 11 | `bda38f5` | Try generalising rkf45               | 10.72 +/- 0.20 ms | 3.49x       | 1.28x       |

Four commits carry the entire 3.49x improvement; everything else is neutral within
noise (~2sigma ~= +/-2%).

## Per-commit analysis

### 02-05 `e861f2f`, `5a5c677`, `011e699`, `cd9998a`  -  formatting/clippy: neutral (as expected)
Pure style changes: rustfmt, `.gitignore`, clippy lints (needless returns, references),
comment format. No change in generated code paths; timings identical within noise.

### 06 `d9f262e`  -  Migrate to num-dual: **2.23x faster**
The big one. Upstream computes the Newton Jacobian by **complex-step
differentiation**: the entire `State` is `Complex64`-valued, so *every* arithmetic
operation in the ODE pays complex cost even when only the real part matters, and
transcendental functions are the killer  -  Sutherland viscosity's `powf`/`sqrt` on a
complex number go through polar form (`atan2`, `exp`, `ln`, `sin`, `cos`), versus a
single real `sqrt`/`powf` plus one multiply-add for a dual number. Switching
`State<Complex64>` -> `State<Dual64>` keeps the same two-integrations-per-Newton-iteration
structure but replaces complex arithmetic with the much cheaper dual-number chain rule.

### 07 `630ec6e`  -  clippy: neutral
Nine-line lint cleanup, no codegen effect.

### 08 `f0bf71e`  -  first_derivative via nested duals: **0.94x (6% regression)**
Deleted the hand-derived analytic `density_viscosity_product_derivative` (and
`soft_max_derivative`) in favour of `num_dual::first_derivative(|g| ..., g)` inside the
ODE. That evaluates the property function on `Dual<T>`  -  dual numbers *nested over the
already-dual state type*  -  so each call does roughly double the arithmetic of the
hand-written derivative formula. A deliberate correctness/maintainability trade
(the analytic derivative can silently drift from the function it differentiates);
the cost was ~1 ms here.

### 09 `30dbb0b`  -  fmt tweaking: neutral
rustfmt config only.

### 10 `0286b68`  -  Both derivatives in one pass: **1.30x faster**
Replaced two separate `Dual64` integrations per Newton iteration (one for d/dfdd, one
for d/dgd) with a single `DualSVec64<2>` integration carrying both derivative
components at once. The function value (`re`) and all the control flow, `soft_max`,
`sqrt`, `powf` real-part work is computed **once instead of twice**; only the cheap
linear derivative-propagation work still scales with the number of eps components.
Theoretical ceiling ~2x on the Jacobian passes; realised 1.30x overall because the
wider 3-lane state (144 bytes) costs more per step and the final plain-`f64`
integrations don't benefit.

### 11 `bda38f5`  -  Generalised rkf45: **1.28x faster**
Three compounding effects in the stepper:
1. **`fn` pointer -> `impl Fn`**: the ODE callback is now monomorphized and statically
   dispatched, guaranteeing the (large, dual-heavy) `self_similar_ode` inlines into all
   six `k`-stage evaluations instead of relying on LLVM devirtualizing a function pointer.
2. **Division -> multiplication**: `3.0 * h * k1 / 32.0` became `k1 * h * (3.0 / 32.0)`  - 
   the tableau divisions on every State component (18 lanes each) became compile-time
   constants multiplied through, and FP division is ~4x the latency of multiplication.
3. Minor: the reassociated arithmetic also changed convergence, 18 -> 16 Newton
   iterations (identical converged results), which is itself ~11% fewer integrations.

Note: this commit also dropped `.abs()` from the returned error estimate (the generic
`S: Add + Mul<f64>` bound can't express it). Harmless today since the error is
discarded, but it matters if adaptive stepping is added.

### 12 working tree  -  scan -> for loop: neutral in this run
Earlier measurements in isolation showed the `once().chain(scan()).collect()` shape
costing ~6% versus a plain loop (external iteration through `Chain`/`Scan` with state
behind `&mut`); in this interleaved run the difference is below noise (10.73 vs
10.72 ms). Either way the for-loop version is never slower, and it removes the
prepend-first-value wart.

## Takeaways
- The 3.49x total = cheaper derivative arithmetic (complex -> dual, 2.23x) x fewer
  passes (batched duals, 1.30x) x a tighter inner loop (inlined + division-free
  stepper, 1.28x), with a 6% correctness tax from nested-dual property derivatives.
- Remaining ideas: recover the `f0bf71e` tax by caching `density_viscosity_product`
  evaluations, and restore `|err|` handling before moving to adaptive stepping.
