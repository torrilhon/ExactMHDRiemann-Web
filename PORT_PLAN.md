# Port Plan: ExactMHDRiemannSolver.jl → TypeScript web tool

Snapshot port of the Julia package to TypeScript. It runs entirely in the browser and is hosted on GitHub Pages.

| | |
| --- | --- |
| Source | [torrilhon/ExactMHDRiemannSolver](https://github.com/torrilhon/ExactMHDRiemannSolver), commit `6da1c024b8e5af6975dd66569e38b936424cb0e9` (2026-10-01), version 0.1.0. Tag this commit in the Julia repo (e.g. `v0.1.0-web-snapshot`). |
| Numerics pin | The Julia repo has no `Manifest.toml` (it is git-ignored) and `[compat]` only fixes major versions (`NonlinearSolve = "4"`, `OrdinaryDiffEqVerner = "2"`, `Roots = "3"`, …). The tag alone therefore does **not** fix the numerics. The real pin is `reference/Manifest.toml`, committed in the web repo (§5.1). |
| Target | Existing repository [torrilhon/ExactMHDRiemann-Web](https://github.com/torrilhon/ExactMHDRiemann-Web) (MIT license already in place), static site on GitHub Pages |
| Policy | One-time snapshot; the two code bases may diverge later. The Julia code is the reference for correctness of this port only. |
| Goal | Same robustness and accuracy as the Julia original: same return codes, same refusals, physical results to ~1e-9. |
| First milestone | The example of the 2002 report (`examples/paper.toml`) solved end to end in the browser. It contains every wave kind of the regular solver (fast fan, rotation, slow fan, contact, slow shock, rotation, fast shock) and so exercises the whole regular pipeline at once. |

---

## Status

| Phase | State |
| --- | --- |
| 0. Setup | done: scaffold, CI, Pages, `reference/` pinned to `6da1c02` (`reference/BASELINE.md`) |
| 1. Golden data | done: 11 golden files, 8021 reference solves (`test/golden/README.md`) |
| 2. Numerical toolbox | done: checked math with exact fma, Julia formatting, linalg, Brent (bit-identical to Roots.jl), Vern9 engine with dense output and step replay, Moré trust region, branch-aware FD Jacobian (83 unit tests) |

Findings so far (details in `reference/BASELINE.md`):

- **Julia bug, fixed on a branch.** Two valid inputs just above the Bt thresholds return `CheckFailed`: the checker reports dense-output noise of a vanishingly weak fan as a mismatch. Fix `0ca567a` on branch `claude/elegant-edison-fb6a3w` of the Julia repo changes only these two results (verified on all 8021 solves). Pending: merge and re-pin the snapshot. The TS port includes the fixed checker; until the re-pin, the two cases are listed in `DIFFERENCES.md`.
- **Time limit.** Five Julia homotopy solves need 6.8–8.3 s, above the default 5 s, so the golden data (time limit off) are not what Julia returns with default options for these cases. Confirms §4.2, item 4.
- **Domain errors** in the residual appear only for |Ψ| ≳ 1e4; the `BIG` path is exercised by 130 golden points.
- **Julia uses FMA.** `evalpoly` is Horner with `muladd`, a hardware fma on the reference machine; the port emulates fma exactly.
- **The trust region is Moré's, not dogleg** (§4.2).
- **Reconstructed scans** match the README, except a > c_A small-Bt settings at 1e-4 with 28/30 instead of ≥ 29/30.

---

## 1. Scope of the first version

In scope:

- All problem classes the Julia snapshot solves: the regular solver, and the quasi-Euler solver for Bn = 0 (any Bt, including Bt = 0). The regular solver comes first (§6); the quasi-Euler path is added afterwards, on top of the finished infrastructure.
- A form for the left and right states (ρ, vx, vy, vz, Bx, By, Bz, p), γ, and the output parameters t, x-range and n.
  - An "advanced" panel with the `SolverOptions` that matter for experiments: at least `bn_euler` (validation set 11 needs `bn_euler = 0`), `time_limit`, `homotopy`, and an optional start vector `guess` (5 path variables, canonical frame).
- Presets: the report example (`examples/paper.toml`), Brio–Wu (`examples/briowu.toml`) and Sod's problem with Bn = 0 (quasi-Euler path, Toro's test 1).
- Load a problem file (`.toml`, the format of `problem_load`) into the form. Use a small TOML parser such as `smol-toml`, with the same defaults and the same checks as `problem_load` (8 numbers per state, `output.x = [xmin, xmax]`). Error messages name the uploaded file, as Julia's name the path.
- A **Compute** button. The solve runs in a Web Worker so the page stays responsive.
- A result panel:
  - `retcode`, `reason` and `method`
  - residual and the worst check error (`check.maxerr`), plus check messages if any
  - the wave table: kind, left and right speeds, and the state right of each wave
  - static preview plots of the 8 primitive variables over x at time t
- Downloads:
  - `solution.csv` with Julia's `write_csv` header and columns `x,rho,vx,vy,vz,Bx,By,Bz,p`, numerically equal to Julia's output (see §4.7 for number formatting)
  - `solution.json` with the input, options, retcode, Ψ, wave table, check report, snapshot commit and port version
  - `problem.toml`, the input in the format of `problem_load`, so any result can be reproduced in Julia with `problem_load` + `solve` or `bin/riemann.jl`

Out of scope for now: live sliders, animation, problem gallery, intermediate/compound waves (not in the Julia v0.1 either).

---

## 2. What has to be ported

The Julia package has 1,287 lines of source in 12 files under `src/` (including `perpendicular.jl`, `problemfile.jl` and the module file). The numerical "tools" it relies on come from external packages, and **these are the real porting work**. The MHD-specific code is mostly plain arithmetic and translates almost line by line.

| Julia dependency | Where it is used | TypeScript replacement |
| --- | --- | --- |
| `NonlinearSolve.TrustRegion` + `AutoForwardDiff` | `nl_solve` (5×5 system, 60 iterations, `abstol = 1e-12`) | Own trust-region (dogleg) solver, ~200 lines (§4.2) |
| `ForwardDiff` (Jacobian, `ift_polish`, `ForwardDiff.derivative`) | Jacobian of the residual; derivatives of roots via the implicit function theorem | Finite-difference Jacobian with branch-aware stencils and frozen ODE step sequences; optional forward-mode duals later (§4.3) |
| `OrdinaryDiffEqVerner.Vern9` (dense, tol 1e-12, `maxiters = 100_000`) | Fast and slow fan ODEs (`_fan_solve`); quasi-Euler fan (`perp_fan`, 1 equation; OrdinaryDiffEq's default `maxiters` is also 100,000) | Own explicit RK engine with Vern9 tableau and dense output (§4.1) |
| `OrdinaryDiffEqVerner.Vern7` (tol 1e-11) | Independent eigenvector integration in `check.jl` | The same Vern9 engine at tol 1e-11, end value only. No second tableau (§4.1). |
| `Roots.Brent` (`find_zero`, `xatol = 0`, `xrtol = 4eps()`; `xatol = 1e-14` in `sample`) | `largest_root_in`, `slow_limit`, `sample`; in `perpendicular.jl`: fan end density on `(1e-300, ρ₀)`, shock compression on `(1, κ(1 − 1e-14))`, P* in the outer bracket | Brent with Roots.jl's stopping tolerances (§4.4) |
| `LinearAlgebra.eigen` (nonsymmetric 7×7), `opnorm(A, Inf)`, `I` | `fan_eigen_rhs`, `fan_pointwise_defect` | `ml-matrix` `EigenvalueDecomposition`, or inverse iteration (§4.5) |
| `StaticArrays.SVector/SMatrix` | everywhere | Plain `number` fields, small tuples or `Float64Array` (§3) |
| `evalpoly`, `hypot`, `atan(y, x)`, `tanh`, `exp`, `log`, `sqrt` | kernels | `Math.*` with checked wrappers (§4.6) |
| `TOML`, `Printf`, `Logging` | `problem_load`, CLI, display, silencing | `smol-toml` (parse and write) for problem files; a Julia-compatible number formatter (§4.7); no logging |

---

## 3. Target architecture

### 3.1 Repository layout

```
ExactMHDRiemann-Web/
  src/
    solver/                     # pure numerics, no DOM; runs in Node and in a Worker
      types.ts                  # RetCode, SolverOptions, HState, Ctx, Wave, FanData, Frame, RiemannProblem, RiemannSolution
      math.ts                   # DomainError, sqrtD, logD, powD, hypot, evalpoly, sign helpers
      linalg.ts                 # 5x5 LU with partial pivoting, mat-vec, inf-norms
      brent.ts                  # Brent root-finder (Roots.jl tolerances)
      rk/
        engine.ts               # adaptive explicit RK, error control, dense output, frozen-step replay
        vern9.ts                # tableau + interpolation coefficients (generated, see §4.1)
      trustregion.ts            # dogleg trust-region for small dense systems
      jacobian.ts               # finite-difference Jacobian (and later dual numbers)
      eos.ts  canonical.ts  shocks.ts  fans.ts  side.ts  check.ts  perpendicular.ts  solve.ts  output.ts  problemfile.ts
      format.ts                 # Julia-compatible repr of Float64, Julia-compatible range(a, b; length)
      index.ts                  # public API: solve, sample, wavetable, toCSV, toJSON
    worker/solver.worker.ts     # message protocol: {problem, options, output} → result
    ui/                         # form, presets, result panel, plots, downloads
  test/
    golden/                     # JSON files generated once by Julia (committed)
    unit/                       # kernel tests ported from runtests.jl
    golden.*.test.ts            # comparison against Julia
  reference/                    # Julia project that produces the golden data
    Project.toml  Manifest.toml # Manifest committed: this is the numerics pin
    generate_golden.jl          # driver; parts in golden/*.jl
    golden/scans.jl             # new: README's small-Bt and 3,800-problem scans (not in the Julia repo)
    BASELINE.md                 # recorded results of the reference (Phase 0/1)
  index.html  vite.config.ts  tsconfig.json  package.json
  .github/workflows/ci.yml      # test + build + deploy to Pages
```

### 3.2 Tooling

- **TypeScript** in `strict` mode with `noUncheckedIndexedAccess`, so array accesses that may be out of bounds are compile errors.
- **Vite** for the dev server and production build, **Vitest** for tests (both run in Node).
- **Plots:** uPlot (small and fast) or Plotly (heavier, has export built in). For static previews, uPlot is enough.
- **CI** (GitHub Actions): `npm ci` → `tsc --noEmit` → `vitest run` (unit + golden) → `vite build` → deploy `dist/` to Pages.
  - Julia is not needed in CI: the golden files are generated once and committed.
  - If the repository is private, GitHub Pages needs a paid plan; check before Phase 0.

### 3.3 Mapping conventions

- **`HState{T}`** → `interface HState { rho; u; p; bt; phi; vt0; vt1 }` with all fields `number`. Flat fields are faster than nested arrays in JS engines.
- **`SVector{5}`** (Ψ, residual) → `Float64Array(5)` or `[number, …]` tuples. **`SVector{8}`** primitive states → `Float64Array(8)`, with index constants `RHO = 0 … P = 7` to avoid off-by-one errors from Julia's 1-based indexing.
- **Symbols** (`:fast_shock`, `:tangential`, `:quasi_euler`, …) → string-literal union types. `RetCode` → `"Success" | "InvalidInput" | "Unsupported" | "RegularLimit" | "NoConvergence" | "CheckFailed"`.
- **`σ::Int`** → `number`, always ±1.
- **Multiple dispatch:**
  - `speeds(h, ctx)` / `speeds(W, γ)` / `speeds(ρ, p, bt, Bn, γ)` → `speeds`, `speedsH`, `speedsW`
  - `sample(sol, ξ)` / `sample(sol, x, t)` → `sampleXi`, `sampleXT`
- **`try … catch err; err isa DomainError || rethrow()`** → `catch (e) { if (!(e instanceof DomainError)) throw e; … }`.
  - This must be kept exactly. Programming errors (TypeError, RangeError) must surface, never be turned into "no convergence".
- **`fvalue(...)`** (strip ForwardDiff partials) → identity. With a finite-difference Jacobian every quantity is differentiated, including the few that Julia deliberately cuts (`yhi`, returned speeds); this is harmless because none of them changes the residual's value.
- **Keep file boundaries and function names** of the Julia code. A reviewer can then compare `shocks.jl` and `shocks.ts` side by side. Add a comment header to every TS function naming the Julia function it mirrors.

---

## 4. The numerical tools in detail

### 4.1 ODE integration (fans and checker)

**What the Julia code does**

- `_fan_solve` integrates the fast fan (3 equations) or slow fan (4 equations) on τ ∈ [0, 1]:
  - Vern9, `abstol = reltol = 1e-12`, `maxiters = 100_000`
  - with `record = true`: dense output (`save_everystep = true`)
  - an unsuccessful retcode becomes `DomainError`
- Dense output is used by:
  - `fan_state(f, τ)`, called inside Brent by `sample` (one root find per output point inside a fan)
  - `fan_pointwise_defect`, a 4th-order finite difference with h = 1e-3 of the dense solution, accepted at `100·check_tol = 1e-6`
- `check.jl` integrates the eigenvector field with Vern7, tol 1e-11, from ρ₀ to ρ₁ (end value only).

**Required accuracy of the dense output: ~1e-10, not 1e-12.** In `fan_pointwise_defect` an interpolation error ε is amplified by about 1/h = 1e3, so ε = 1e-10 gives a defect of ~1e-7, well below 1e-6. `sample` is compared at 1e-10 (§5.3). Any high-order interpolant at tol 1e-12 meets this.

**Plan**

1. **Write a generic adaptive explicit RK engine.**
   - Input: tableau `(c, A, b, b̂)`, optional interpolation matrix for extra stages.
   - Error norm: RMS over components of `err_i / (abstol + reltol·max(|y_i|, |y_new_i|))`, the OrdinaryDiffEq default internal norm (see the [DiffEqBase API docs](https://sciml.github.io/OrdinaryDiffEq.jl/stable/api/diffeqbase/)).
   - Step-size controller: a standard controller (PI or plain I with safety factor, `qmin`/`qmax`) and a standard initial-dt heuristic. **No attempt to replicate OrdinaryDiffEq's controller exactly**: at tol 1e-12 any correct controller agrees with Julia far below the acceptance limits of §5.3.
   - **Frozen-step replay:** besides the adaptive mode, the engine can re-integrate with a given step sequence (the one recorded by a previous adaptive run). This is what the Jacobian uses (§4.3).
   - `maxiters` and non-finite stage values → `DomainError`, as in Julia.
2. **Vern9 coefficients and the 9th-order interpolant**
   - Primary source: `lib/OrdinaryDiffEqVerner/src/verner_tableaus.jl` in [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl/tree/master/lib/OrdinaryDiffEqVerner/src), MIT license. The interpolant lives in `interpolants.jl` / `verner_addsteps.jl` (extra stages used only for dense output).
   - Secondary source, to cross-check digits: [Jim Verner's page](https://www.sfu.ca/~jverner/), e.g. the ["most robust" 9(8) pair](https://www.sfu.ca/~jverner/RKV98.IIa.Robust.000000351.081209.FLOAT6040OnWeb). Make sure the variant matches what OrdinaryDiffEq calls Vern9 ("efficient" vs "robust").
   - Generate `vern9.ts` with a small script, never by hand. Copy coefficients as decimal strings with ≥ 17 significant digits.
3. **The checker uses the same engine.** `fan_eigen_integrate` runs Vern9 at tol 1e-11 without dense output instead of Vern7. The checker is an independent consistency check at `check_tol = 1e-8`; its integrator does not need to match Julia's, only to be accurate. This saves a second tableau and its tests.
4. **Verify the engine on its own before using it in the solver:**
   - Order tests: on y' = −y and y' = cos(t)·y, the global error should scale with the expected order as tol is lowered.
   - Interpolant test: the dense solution at 1,000 random τ agrees with a re-integration to that τ, to ≲ 10·tol.
   - Tableau consistency: row sums `Σ_j a_ij = c_i`, `Σ b_i = 1`, order conditions up to order 3 as a cheap typo detector.
   - Frozen-step replay reproduces the adaptive run bit for bit when given its own step sequence.
   - Compare with Julia: golden end values and dense values of `fast_fan` / `slow_fan` for ~200 random states (§5).
5. **Fallback, possibly the simpler first choice:** Hairer's DOP853 with its 7th-order dense output.
   - Sources: [Hairer's codes](https://www.unige.ch/~hairer/software.html), a modern, readable refactoring in [jacobwilliams/dop853](https://github.com/jacobwilliams/dop853), and Hairer–Nørsett–Wanner, Solving ODEs I, Sec. II.6.
   - At tol 1e-12 its dense output easily meets the ~1e-10 requirement above. Decide in Phase 2 by whichever tableau and interpolant is quicker to get verifiably right.
   - Avoid generic npm ODE packages (e.g. RK45 / Cash–Karp ports): wrong order for 1e-12, and usually no high-order dense output.

### 4.2 Nonlinear solver (5×5)

**What the Julia code does**

- `nl_solve` calls `NonlinearSolve.TrustRegion(autodiff = AutoForwardDiff())` on `safe_residual`, with `abstol = 1e-12` and `maxiters = 60`.
- `safe_residual` maps any `DomainError` or non-finite value to the vector `(1e6, …)`, so the trust region simply rejects such steps.
- Solver exceptions are swallowed, except a list of programming errors.
- Acceptance is decided separately, by `accept(r) = r ≤ max(1e3·abstol, 1e-10)` on the ∞-norm. **Bit-for-bit agreement with NonlinearSolve's internal stopping is not required.**
- Around it: a fixed list of start vectors (`starts(α)`), an optional user guess, then the homotopy in λ:
  - predictor–corrector with step halving and doubling (`h` from 0.05, ≤ 0.25, stop below 1e-4)
  - `homotopy_maxsteps = 400`
  - wall-clock `time_limit = 5 s`

**What `TrustRegion()` actually is at the pinned version** (read from the sources in Phase 2; the NonlinearSolve docs describe an older default):

- descent `MoreTrustRegionDescent` (not dogleg): the subproblem `min ‖J p + F‖, ‖p‖ ≤ Δ` is solved nearly exactly by Moré's iteration on the damping λ (MINPACK `lmpar`), D = I, θ = 1e-4, at most 10 λ updates, λ warm-started; Gauss-Newton step when it is inside the region and J has full numerical rank (pivoted QR)
- radius update `RadiusUpdateSchemes.More` (MINPACK `lmder`): accept at ρ > 1e-3; ρ < 1/4 → Δ = ¼ min(Δ, 10‖p‖); ρ ≥ 3/4 or λ = 0 → Δ = 2‖p‖; initial Δ = ‖u₀‖₂; predicted reduction ½‖Jp‖² + λ‖p‖²
- Jacobian recomputed only after accepted steps; `max_shrink_times = 32`
- termination `AbsNormSafeBestTerminationMode(‖·‖∞, max_stalled_steps = 32)`: the best iterate is returned; patience (100 steps) never applies with `maxiters = 60`

On 16 classic test problems (Moré–Garbow–Hillstrom) the TS port (`trustregion.ts`) returns the same retcodes as Julia and, in 13 cases, exactly the same number of steps (`test/unit/trustregion.test.ts`).

**Plan**

1. ~~Implement a dogleg trust-region~~ Done as a port of the Moré trust region above (`trustregion.ts`). Original plan, for reference:
   - merit ½‖F‖²
   - Gauss–Newton step from LU with partial pivoting (`linalg.ts`)
   - Cauchy step, dogleg path
   - ratio ρ = actual/predicted reduction
   - radius update with exactly the constants listed above
   - termination on `‖F‖∞ ≤ abstol`, `maxiters`, or `max_shrink_times` consecutive shrinks, plus the stall rules of the termination mode found above
   - return the iterate selected by that termination mode (likely the best one seen)
   - singular J: fall back to the Cauchy (steepest-descent) step instead of throwing
2. Reference algorithm for robustness details (scaling, singular Jacobians): MINPACK `hybrj` (Powell hybrid). See the [netlib source](https://www.netlib.org/minpack/) and [documentation](https://www.math.utah.edu/software/minpack/minpack/hybrj.html). Useful for reading, not for porting wholesale.
3. Port `starts`, `homotopy`, `accept`, `limit_reason` and `_solve` **literally**. They are plain control flow and carry most of the robustness. Use `performance.now()` for the time limit. The Worker gets its own extra hard timeout (e.g. 20 s) as a last safety net.
4. **The time limit is a performance-dependent decision.** Without ForwardDiff's single dual pass, every Jacobian costs 10 residual evaluations, each with up to four Vern9 fan integrations and possibly `slow_limit`'s 64-point scan. The homotopy may run 5–10× slower than in Julia, and a 5 s wall-clock limit would then change retcodes.
   - For all golden comparisons, run with `time_limit = Infinity`, so the algorithms are compared, not the machines.
   - Measure the homotopy times of the TS port separately (stress and small-Bt sets), then set the shipped default so no case that succeeds without a limit is lost. Optionally add a step-count limit, which is machine-independent.
5. **Robustness benchmark, not just accuracy:** the stress and small-Bt sets (§5) must reach the Julia success rates. If the TS trust region does worse on some cases, compare iteration histories on those cases first. Don't change tolerances.

### 4.3 Derivatives: replacing ForwardDiff

ForwardDiff is used in two ways:

- (a) the full 5×5 Jacobian of the residual for the trust region
- (b) `ift_polish`, which makes derivatives of root-finder outputs (cubic root in `fast_shock`, `slow_limit`) exact by the implicit function theorem. For Float64 inputs it is simply **one Newton polish step**, which also improves the value.

JavaScript has no operator overloading, so the Julia trick "same code, dual numbers" does not carry over directly.

**Phase 1: finite differences (recommended start)**

- **Central differences** per component: `h_j = ε^{1/3} · max(1, |Ψ_j|)` ≈ 6e-6·max(1, |Ψ_j|), i.e. 10 residual evaluations per Jacobian. Truncation error is O(h²) ≈ 1e-11. `runtests.jl` already checks that ForwardDiff and central differences with h = 1e-6 agree to 1e-6 on the paper example, so the approach is known to work there.
- **Frozen ODE step sequences.** An adaptive integration is a smooth function of its parameters only while the controller makes the same accept/reject decisions and the same last-step clipping. Jumps of order ode_tol (1e-12) divided by h would give relative Jacobian errors up to ~2e-7. To avoid this, the base residual evaluation records the step sequence of every fan integration, and the perturbed evaluations replay those sequences with the frozen-step mode of the engine (§4.1). This is internal numerical differentiation (Bock): the FD Jacobian is then a smooth derivative of one fixed discrete map. It is cheap and should make Phase 2 unnecessary.
- **Branch-aware stencils.** `build_side` switches parametrization at branch points:
  - fast wave at ψf = 0: shock (`ψ = √D − √D_f`) for ψf > 0, fan (`ℓ = log(Bt/Bt₀)`) for ψf ≤ 0
  - slow wave at ψs = `SLOW_EPS` = 1e-9: shock (B̂t drop Δ through tanh) above, fan (arc length ς in (s, B̂t)) below

  The wave curve is continuous across the switch, but its parametrization speed is not, so the residual has a kink there. ForwardDiff always differentiates exactly one branch. A central stencil with |ψ − switch| < h averages the two slopes. This happens for vanishing waves (the "single-wave problems" testset) and near the start vectors (ψ = ±0.01 is far enough, but iterates can move closer). **Rule: if the central stencil would cross a switch, use the one-sided difference on the base point's side.**
- **Domain edges:** if a perturbed point returns the `BIG` sentinel, use the one-sided difference on the other side. If both sides fail, mark the column as unreliable and let the trust region shrink. Without this, Jacobians near the slow-shock limit or vacuum blow up.
- `ift_polish` stays as a value-level Newton step, `x0 − f(x0)/f'(x0)`, with `f'` evaluated by an analytic derivative (the cubic: `evalpoly` of the derivative coefficients) or a central difference (`slow_defect`).
- Optional improvement: Broyden updates between full FD Jacobians to cut residual evaluations. Measure first; the current cost is likely fine.

**Phase 2 (only if Phase 1 misses the robustness targets): forward-mode duals**

- A small `Dual` class with explicit functions (`dadd`, `dmul`, `dsqrt`, …) and kernels written against a tiny numeric interface.
- Cost: the kernels in `eos/shocks/fans/side` must be written in that style, and the RK engine must run on duals. That doubles the review effort.
- Alternative: hand-derived Jacobian blocks for the shock kernels (closed-form polynomials), with FD only through the fans.
- Decide based on the comparison of iteration counts and success rates against Julia.

### 4.4 Scalar root finding (Brent)

- Implement Brent with **Roots.jl's stopping tolerances**, using Roots.jl's `Brent` (MIT) as the reference for its termination tests. The call sites rely on:
  - `xatol = 0.0, xrtol = 4eps()` → iterate until the bracket is within ~4 ULP of the root (`largest_root_in`, `slow_limit`, all three calls in `perpendicular.jl`)
  - `xatol = 1e-14` in `sample`
  - an error when the bracket has no sign change
  - early return when an endpoint is an exact zero (check Roots.jl's handling of `f(a) == 0` and `f(b) == 0`)
- **Not required:** identical iteration sequences or iteration counts. With brackets this tight, any correct Brent with the same tolerances agrees with Roots.jl to a few ULP.
- Generic JS versions (e.g. the [Wikipedia-derived gist](https://gist.github.com/ryanspradlin/18c1010b7dd2d875284933d018c5c908) or [uniroot.js](https://gist.github.com/borgar/3317728)) are fine for cross-checking but have different defaults.
- Tests:
  - polynomials with known roots, to a few ULP
  - the cubic `largest_root_in` on golden coefficient sets from Julia, including near-switch-on cases (A → 0), near-double roots, and Bn/√p down to 1e-10 (vanishing leading coefficient, spurious far root)
  - very wide brackets as in `perp_fan`, `(1e-300, ρ₀)`: Brent must reach the root in a reasonable number of steps (bounded, e.g. ≤ 200)

### 4.5 Linear algebra

- **5×5 LU** with partial pivoting and **7×7 mat-vec**: hand-written, a few dozen lines, unit-tested against hand-computed cases.
- **`opnorm(A, Inf)`**: maximum absolute row sum.
- **Nonsymmetric 7×7 eigenproblem** (`fan_eigen_rhs` in the checker, called at every RK stage of the independent fan integration), two options:
  - A. [`ml-matrix`](https://github.com/mljs/matrix) `EigenvalueDecomposition` (JAMA port of EISPACK `hqr2`, handles real nonsymmetric matrices; MIT). Closest to Julia's LAPACK `geev`. The code normalizes the eigenvector by `r[1]`, so sign and scale conventions don't matter.
  - B. Inverse iteration on `(A − λI)` with λ = u ± c_f / u ± c_s from `speeds`. Cheaper and dependency-free, but it ties the checker's eigenvector more closely to `speeds`. The check is already allowed to share `speeds`, so this is acceptable.
  - Recommendation: A for faithfulness. Keep B as a test-only cross-check.
  - Note: `sortperm(real.(E.values))[fam]` picks the fam-th smallest eigenvalue (sorted: u−cf, u−cA, u−cs, u, u+cs, u+cA, u+cf; fam ∈ {1, 3, 5, 7}). Copy this sort order exactly, with a numeric comparator.
  - With Bn ≈ 0 (quasi-Euler fans) the eigenvalue u is fivefold. Fam 1 and 7 are still simple, but test the chosen eigen solver on such matrices.

### 4.6 Floating-point semantics: the subtle part

These differences between Julia and JavaScript would silently change behaviour:

| Topic | Julia | JavaScript | Action |
| --- | --- | --- | --- |
| `sqrt(x)`, x < 0 | throws `DomainError` | returns `NaN` | `sqrtD(x)` that throws `DomainError` (see the rule below) |
| `log(x)`, x < 0 | throws `DomainError` | `NaN` | `logD` |
| `x^y`, x < 0, y non-integer | throws `DomainError` | `NaN` | `powD` for `(ρ/ρ₀)^γ` in `perpendicular.jl` |
| `x^2`, `x^3` (literal) | `x*x`, `x*x*x` (`literal_pow`) | `x ** 3` may call `pow` | write `x*x`, `x*x*x` to get the same rounding |
| `evalpoly(x, (c0, c1, c2, c3))` | Horner with `muladd`, coefficients ascending; **a hardware FMA** on x86-64/ARM64 | no fma | `evalpoly` with an exact fma emulation (`math.ts`), bit-identical to Julia |
| `atan(y, x)` | two-argument atan | `Math.atan2(y, x)` | wrapper `atan2` |
| `hypot` | robust | `Math.hypot` robust, slower, last-bit differences possible | fine; inline only if profiling says so |
| `sort(v)` | numeric | **lexicographic by default** | always `sort((a, b) => a - b)`; matters in `largest_root_in` (`sort([lo; crit; hi])`) and for `sortperm` in `fan_eigen_rhs` |
| `max(NaN, x)` | `NaN` | `Math.max` → `NaN` | same, but guard any `>`/`<` comparisons that rely on NaN behaviour |
| `eps()` | 2.22e-16 | `Number.EPSILON` | constant |
| Division by 0, overflow | IEEE Inf/NaN | IEEE Inf/NaN | same |
| `Math.exp/log/tanh/sin/cos` | openlibm, ≤ 1 ULP | engine-specific, usually ≤ 1 ULP, **can differ between browsers** | compare with tolerances, never bitwise; run the golden tests in Node **and** in at least one browser (Playwright) |
| 1-based indexing | `W[1] … W[8]` | `W[0] … W[7]` | named index constants; one dedicated review pass for indices |

**Rule: every `sqrt` and `log` in the solver becomes `sqrtD` / `logD`, without exception.** No judgment about which arguments "could" be negative; the cost is nil. Every `throw(DomainError(...))` in the Julia code becomes a `DomainError` in TS. A grep in CI (no bare `Math.sqrt` / `Math.log` / `Math.pow` under `src/solver/`) keeps this enforced.

Why it matters even though `safe_residual` also catches non-finite values: the residual itself would be safe with `NaN`, but the control flow around it would not. The decisive case is `slow_limit`:

- In Julia, a `DomainError` in the first evaluation `fv(Δprev)` is not caught there. It escapes to `safe_residual`, which returns the `BIG` vector.
- In JS with a bare `NaN`, `fprev >= 0` is false and the scan silently continues, returning some Δmax from a meaningless branch.

### 4.7 Output formats: CSV and JSON

**Number formatting.** Julia's `repr(Float64)` and JS `String(number)` are both shortest round-trip representations, but they print differently:

| value | Julia `repr` | JS `String` |
| --- | --- | --- |
| 1 | `1.0` | `1` |
| 1e-5 | `1.0e-5` | `0.00001` |
| −0 | `-0.0` | `0` |
| 1e21 | `1.0e21` | `1e+21` |

`format.ts` provides `juliaRepr(x)` with Julia's rules (always a decimal point; scientific notation outside Julia's fixed-notation range, with `e-5` style exponents; `-0.0`), tested against a golden list of values printed by Julia. CSV files then agree character for character when the numbers agree.

**The x grid.** Julia's `range(a, b; length = n)` is a `StepRangeLen` with TwicePrecision arithmetic, so its points can differ by 1 ULP from a naive `a + i*(b − a)/(n − 1)`. `format.ts` provides `juliaRange(a, b, n)`, which reproduces Julia's values, tested against Julia's output for the grids of the examples and some awkward ones (e.g. `[-0.5, 0.5]`, `[-1, 3]` with n = 7, 2001, 1000).

**JSON.** JSON has no `Inf` or `NaN` (Julia stores `residual = Inf` on failure). Encode non-finite numbers as the strings `"Inf"`, `"-Inf"`, `"NaN"`, both in `solution.json` and in the golden files, and document this in the JSON schema.

### 4.8 Points of attention in the source

Details of the snapshot that are easy to get wrong in a port:

- **`validate_input` order.** Checks run in this order: non-finite → γ → density → pressure → `:bn_jump` → `:extreme_ratio` → Bn/√p ≤ `bn_euler` on both sides returns `(Success, :quasi_euler)` and skips the Bt thresholds → `:switch_on_off`. The order decides which reason is reported when several apply. Keep it literally.
- **Quasi-Euler dispatch.** `_solve` calls `solve_perpendicular` for `:quasi_euler` and maps a `DomainError` there to `NoConvergence` with `frame = nothing` and `ctx = nothing`. The result type, `wavetable`, and the JSON export must allow a missing frame.
- **Quasi-Euler kernels** (`perpendicular.jl`):
  - wave kind `:tangential` (u and total pressure p + Bt²/2 continuous), fan family `:fast0` with s = s*·τ and Bt ∝ ρ, `method = :quasi_euler`
  - `ctx.Bn` is the user's `L[5]`, which may be tiny but nonzero; `speeds` then includes it. Port literally.
  - P* by Brent on `(1e-12·min(P_L, P_R), hi)`, with `hi` doubled up to 200 times; `g(lo) < 0` means vacuum (`RegularLimit :vacuum`)
  - `perp_shock` throws `DomainError` when P/P₀ ≳ 1e14 (compression indistinguishable from κ in Float64)
  - the residual reported for quasi-Euler is |u*_L − u*_R| over a fast-speed scale, and Ψ = (P*, 0, 0, 0, 0)
  - `solve_perpendicular` runs the **full** `check_waves`: quasi-Euler fans have kind `:fast_fan`, so they go through `fan_pointwise_defect` and, because `speed_gap ≈ 1`, through `fan_eigen_integrate` with the 7×7 eigen-decomposition; shocks go through `lax_region`. The quasi-Euler path therefore needs the complete checker.
- **Fast-shock cubic.** `largest_root_in(c, lo, hi)` searches only the physical range `0 < y ≤ √(A² + 2(1+X))·(1+1e-12)`. As Bn → 0 the leading coefficient ∝ B² vanishes and a spurious root ~1/B² appears far outside that range. There is no early return for `f(a) == 0`, only for `f(b) == 0`.
- **Slow shock p̂** comes from the momentum balance, `1 − X(v̂ − 1) − (B̂t² − A²)/2`, not from the Hugoniot form (which divides by 1 − κv̂ and breaks down as Bn → 0).
- **`fan_speed`** uses c_s only for the `:slow` family and c_f for `:fast` and `:fast0`.
- **`check.jl`.** `speed_gap` is measured relative to c_f; for the slow family it includes the distance to the contact (`min(c_A − c_s, c_s)/c_f`). The `:tangential` branch checks continuity of u and total pressure, and that the speed equals the flow speed.
- **Problem files.** `problem_load` defaults: γ = 5/3, t = 1, x = [−1, 1], n = 2001. It raises `ArgumentError` for a missing state, a state that is not 8 numbers, and `output.x` not of length 2. The website's file loader should give the same messages.

---

## 5. Validation strategy

The port is accepted only when it reproduces the Julia results across thousands of cases. The reference data are produced **once**, by Julia, and committed as JSON.

### 5.1 Generating the golden data

`reference/generate_golden.jl` writes `test/golden/*.json`.

- The Julia project in `reference/` depends on ExactMHDRiemannSolver at the snapshot tag. Its `Manifest.toml` is committed and **is the pin of all numerical dependencies** (NonlinearSolve, OrdinaryDiffEqVerner, Roots, ForwardDiff, …), since the Julia repo itself has none.
- Every golden file records: the Julia commit, the Julia version, and the versions of all packages from the Manifest. A CI test asserts that these match `reference/Manifest.toml` (commit of ExactMHDRiemannSolver and the package versions).
- All values are written with `repr` (round-trip precision); non-finite values as `"Inf"`, `"-Inf"`, `"NaN"` (§4.7).
- Internal functions are accessed through `const E = ExactMHDRiemannSolver`, as `test/runtests.jl` already does.
- The script also records the run time per solve, for the performance comparison.

| Level | File | Content | Size |
| --- | --- | --- | --- |
| L0 basics | `basics.json` | `speeds`, `conserved`, `flux`, `primitive_jacobian`; canonical frame round trip `from_canonical(_to_canonical(W))` | ~500 random states, γ ∈ {1.4, 5/3, 2} |
| L0 format | `format.json` | `repr` of selected Float64 values; `range(a, b; length = n)` for the example grids and awkward ones (§4.7) | small |
| L1 kernels | `kernels.json` | `fast_shock` for ψ ∈ {1e-10, 1e-4, 0.05, 0.5, 3}; `slow_limit`; `slow_shock_state` at fractions {1e-7, 1e-3, 0.3, 0.9, 0.999}·Δmax; `rotation`; `fast_fan(ψ = −0.7)`, `slow_fan(ψ = −0.4)` end states; `largest_root_in` coefficients, interval and root; `fast_shock` with Bn/√p = 1e-9 (as in `runtests.jl`); `perp_shock`, `perp_fan` for P above and below P₀ | ~300 random states, as in `runtests.jl` |
| L1 dense | `fans_dense.json` | `fan_state` at 50 τ values per fan; `fan_speed` | ~100 fans |
| L2 residual | `residual.json` | `residual(Ψ, UL, UR, ctx)` at random Ψ, including Ψ near domain edges (where Julia returns the `BIG` vector) and Ψ near the branch switches (ψf ≈ 0, ψs ≈ `SLOW_EPS`); the ForwardDiff Jacobian at a subset of these points | ~1,000 evaluations, ~100 Jacobians |
| L3 solves | `solve_*.json` | for each problem: `retcode`, `reason`, `method`, Ψ, residual, wave table, `check.maxerr` | see 5.2 |
| L4 output | `csv_*.json` | `sample` on the `write_csv` grid for the examples and 50 random problems | — |

Two validation sets of the README have **no script in the Julia repo** and must be written first; they are in `reference/golden/scans.jl` (sets 5 and 9 below). Run them in Julia and check that they reproduce the README figures before using them as golden data. If they don't, ask the Julia author for the original scripts rather than adjusting the README claims.

### 5.2 Problem sets for L3

1. **Examples:** `paper.toml` (report Tables 1–2, twist angle 1.5; wave sequence fast fan, rotation, slow fan, contact, slow shock, rotation, fast shock) and `briowu.toml`.
2. **Random set:** the 2,000 problems in `test/data/random_problems.csv` (all `Success` in Julia).
3. **C reference:** `test/data/c_reference.csv` (355 problems, states left and right of the contact, 8 decimals). Julia agrees with it to 1.3e-8; the TS port must as well.
4. **Stress set:** the 1,000 problems of `benchmark/stress.jl` (Julia: 97.8 % `Success`, 2.2 % `RegularLimit`, no `NoConvergence`). Export the generated inputs, because the seeded `MersenneTwister` stream cannot be reproduced in JS.
5. **Small-Bt scan:** the "30 random problems per setting" from the README validation section, both regimes c_A ≷ a. New script (§5.1).
6. **Refusals:** one case per `InvalidInput` / `Unsupported` reason:
   - `:nonfinite`, `:gamma`, `:density`, `:pressure`, `:bn_jump`
   - `:switch_on_off`, `:extreme_ratio` (including a case with Bn = 0 **and** an extreme ratio, which must give `:extreme_ratio`)
   - inputs exactly at the thresholds (the `1 − 1e-12` tolerance in `validate_input`, and Bn/√p exactly at `bn_euler`)
7. **Trivial case:** L = R (`method = :trivial`).
8. **Symmetry set:** a few problems under Galilean shift, transverse rotation, B → −B and x → −x (also tested as properties, §5.4).
9. **Quasi-Euler (Bn = 0):**
   - Sod (γ = 1.4): Toro's p* = 0.30313, u* = 0.92745, ρ*L = 0.42632, ρ*R = 0.26557, shock speed 1.75216
   - the 200 random perpendicular problems of `runtests.jl` (`MersenneTwister(5)`, Bt ∈ {0, 1e-8, O(1)}); export the inputs
   - the README's 3,800 perpendicular problems (Bt/√p from 0 to 1e3, ratios up to 1:10⁴). New script (§5.1).
   - the colliding-flow case `vx = ±1e8` (`NoConvergence`, `method = :quasi_euler`)
10. **Switch continuity:** the first 5 random problems with Bn/√p = 2e-10 (regular solver) and 0.9e-10 (quasi-Euler). Outer wave speeds, u* and P* must agree to 1e-8, as in Julia.
11. **Small normal field, regular solver forced** (`bn_euler = 0`): random problems 15, 18, 27 with Bn/√p ∈ {1e-7, 1e-9, 1e-10}, all `Success`. This is the hardest test of `largest_root_in` and the Rayleigh-line p̂.
12. **Vanishing waves:** the "single-wave problems" of `runtests.jl`, where some path variables sit at a branch switch. This set tests the branch-aware FD stencils of §4.3.

All L3 comparisons run with `time_limit = Infinity` in both Julia and TS (§4.2).

### 5.3 Acceptance criteria

The TS solver iterates to the same `abstol = 1e-12` as Julia, not just to the `accept` threshold of 1e-10, so that both converge as far as the problem allows.

| Quantity | Criterion |
| --- | --- |
| `retcode` + `reason` + `method` | identical for **all** problems in sets 1–3, 6, 7, 9–12. For sets 4–5, identical in ≥ 99.5 % of cases. Every mismatch is listed and explained in `test/golden/DIFFERENCES.md`; `Success` vs `RegularLimit` at a limit is acceptable, a different wave structure is not. |
| Hard rule | **no `CheckFailed`**, and no `Success` with check messages (the thresholds of `check_waves`: `check_tol` for jump conditions, `100·check_tol` for fans) |
| Middle states, wave speeds (primary criterion) | `relerr` as in Julia (relative to max(\|a\|, \|b\|, 1)) ≤ 1e-9, or ≤ 1e4 × the larger of the two reported residuals, whichever is larger |
| Ψ (secondary, diagnostic) | `relerr` ≤ 1e-7; the twist angle Ψ[3] compared modulo 2π; cases with a slow shock near saturation (ψs/Δmax > 3, where Ψ is badly determined) are excluded. A Ψ mismatch with matching states is logged, not a failure. |
| `check.maxerr` | ≤ max(1e-8, 10 × Julia's value). Not a pass/fail criterion on its own: fans are accepted up to `100·check_tol`, and Julia's own worst `Success` on the stress set has 1.95e-8 (`reference/BASELINE.md`). |
| C reference | ≤ 1.3e-8 (the Julia figure) |
| Quasi-Euler | P* relative ≤ 1e-13 against Julia (Brent to 4 ULP); states ≤ 1e-12; fan interiors ≤ 1e-10; Sod within 1e-5 of Toro |
| Switch continuity | ≤ 1e-8, as in the Julia test |
| L0/L1 kernels | ≤ 1e-13 for closed-form kernels; ≤ 1e-10 for fan end states (ODE tolerance) |
| L2 Jacobians | FD vs ForwardDiff ≤ 1e-6 × max(1, ‖J‖), as in `runtests.jl` |
| Sampled profiles | ≤ 1e-10 away from discontinuities. Points within 1e-12 of a wave speed are excluded, since they may fall on different sides. |
| CSV | character-identical to Julia's `write_csv` for the examples (§4.7) |
| Speed | single direct solve ≤ 50 ms median in Node, including sampling 2,001 points (Julia: a few ms). Homotopy times are measured and reported, and determine the shipped `time_limit` (§4.2). |

### 5.4 Independent tests (no golden data needed)

Port the property tests of `test/runtests.jl`:

- Rankine–Hugoniot defect of every shock and rotation kernel < 1e-12 over random states and strengths
- fans vs eigenvector integration < 1e-8
- invariance of the solution under the symmetries listed above
- `validate_input` boundaries

Plus:

- RH, Lax and entropy checks via the ported `check_waves` on every golden solve. These are an internal consistency check that does not trust Julia.
- A **fuzz harness** (`npm run fuzz`) that draws random problems in the domain of `random_problems.csv`, solves them in TS and fails on any `CheckFailed`, any uncaught exception, or any `Success` with check messages. Run nightly or on demand, not on every push.
- Browser run: the L0–L3 golden tests also run in headless Chromium and Firefox (Vitest browser mode or Playwright). This catches engine-specific `Math.*` differences.
- Lint check: no bare `Math.sqrt`, `Math.log`, `Math.pow` or `**` with non-literal exponent under `src/solver/` (§4.6).

---

## 6. Work plan

Estimates are focused working days for one person who knows both languages.

The plan goes **straight to the regular MHD solver**, with the report example as the first end-to-end target. Doing the quasi-Euler path first has little benefit: it needs the complete checker (eigen-decomposition, pointwise fan defect, Lax regions), the canonical-frame code, `validate_input` and the output layer, which is most of the infrastructure, while skipping only the trust region and the regular kernels. The report example, on the other hand, exercises every wave kind of the regular solver in one problem and has published reference values (Tables 1–2, 6 digits).

| Phase | Content | Deliverable | Effort |
| --- | --- | --- | --- |
| 0. Setup | Tag the snapshot commit in the Julia repo; Vite/Vitest/TS strict in `ExactMHDRiemann-Web`; CI skeleton with Pages deploy of a placeholder page; `reference/` Julia project with committed `Manifest.toml` | green CI, placeholder page live | 0.5 d |
| 1. Golden data | `generate_golden.jl` for L0–L4 and all problem sets; write the small-Bt and perpendicular scans (`golden/scans.jl`) and check them against the README figures; commit gzipped JSON | `test/golden/` | 2 d |
| 2. Numerical toolbox | `math.ts` (checked functions), `linalg.ts`, `brent.ts`, RK engine with dense output and frozen-step replay + Vern9 (or DOP853) tableau and interpolant (generated), trust region, FD Jacobian; `format.ts` (`juliaRepr`, `juliaRange`); each with its own unit tests (§4) | toolbox passes standalone tests | 4–5 d |
| 3. Kernels | `types`, `eos`, `canonical`, `shocks` (`largest_root_in`, Rayleigh-line p̂), `fans`, `side` | L0, L1 and L2 golden tests green, including FD vs ForwardDiff Jacobians | 2–3 d |
| 4. First end-to-end: report example | `solve.ts` (starts, homotopy, limits, `_solve`, `validate_input`), `check.ts` (incl. eigen choice), `output.ts` (`sample`, `wavetable`, CSV/JSON); minimal page with the form, the `paper.toml` preset, Worker, wave table, one plot and the CSV download | **first live version**: the report example reproduces Tables 1–2 in the browser and its CSV equals Julia's character for character | 3–4 d |
| 5. Regular solver complete | L3/L4 on sets 1–3, 6, 7, 10–12; Brio–Wu preset | golden tests green on these sets | 1–2 d |
| 6. Robustness | Stress and small-Bt sets; homotopy timing and the shipped `time_limit`; analyse every mismatch; `DIFFERENCES.md`; fuzz harness; browser test run | acceptance criteria of §5.3 met for the regular solver | 2–4 d |
| 7. Quasi-Euler | `perpendicular.ts` on top of the finished toolbox, checker and output; `:fast0` fans in `fan_state`/`sample`; Sod preset; set 9 and the quasi-Euler half of set 10 | quasi-Euler golden tests green | 1–1.5 d |
| 8. Web UI | Full form with validation and advanced options, all presets, problem-file upload (`problemfile.ts`), result panel, uPlot previews of all 8 variables, three downloads, plain-language messages per `reason`, citation and "snapshot of commit …" note | live on GitHub Pages | 2–3 d |
| 9. Release | README (scope, accuracy statement with the measured numbers, how to reproduce in Julia via `problem.toml`), `CITATION.cff` pointing to the report and the Julia repo, tag v0.1.0 | v0.1.0 | 0.5 d |

Total: about 18.5–25.5 working days, **4.5–5 weeks**. Phases 2 and 6 carry most of the uncertainty.

**Order of risk reduction:**

1. The RK engine with dense output and frozen-step replay first (Phase 2). If the Vern9 interpolant causes trouble, switching to DOP853 early is cheap.
2. The FD Jacobian is validated against ForwardDiff Jacobians from the golden data (L2) **before** the trust region is used in the driver, including points near branch switches.
3. The report example end to end (Phase 4) proves every component of the regular pipeline: both fan types with dense output, both shock types, rotations, the trust region, the checker and the output, all against published values.

---

## 7. Risks and mitigations

| Risk | Effect | Mitigation |
| --- | --- | --- |
| Vern9 interpolant coefficients mistyped or wrong variant | wrong `sample` inside fans; pointwise fan check fails | generate from source, consistency and order tests, dense golden tests; DOP853 fallback (its accuracy is sufficient, §4.1) |
| FD Jacobian weaker than ForwardDiff near domain edges | lower success rate on stress / small-Bt sets | one-sided differences at the edges, iteration-history comparison, Phase-2 duals or analytic shock Jacobians if needed |
| FD stencil across a branch switch (ψf = 0, ψs = `SLOW_EPS`) | averaged slopes, slow or failed convergence for vanishing waves | branch-aware one-sided stencils (§4.3), set 12, L2 Jacobian tests near the switches |
| ODE step-size decisions differ between base and perturbed evaluation | noisy Jacobian (~1e-7 relative) | frozen-step replay (§4.3) |
| Trust region returns a different iterate than NonlinearSolve (last vs best) | different `NoConvergence` / homotopy behaviour | read the termination mode at the pinned version; return the same kind of iterate (§4.2) |
| Wall-clock `time_limit` hit in TS where Julia finishes | `NoConvergence` that is only a speed artifact | golden comparisons without time limit; shipped limit set from measured TS times; optional step-count limit |
| `NaN` instead of `DomainError` | residual becomes `NaN`, `slow_limit` scans on, wrong limits | checked `sqrtD`/`logD`/`powD` everywhere, lint check, a test that injects out-of-domain Ψ and expects the `BIG` sentinel (L2 golden data contain such points) |
| Index errors (1-based → 0-based) | wrong components, often subtle | named constants, per-module golden tests, side-by-side review |
| Engine-specific `Math.*` rounding | tiny differences, possibly flipping borderline retcodes | tolerances, browser test run, documented differences |
| Number formatting / range differences | CSV not identical to Julia's | `juliaRepr`, `juliaRange`, golden format tests (§4.7) |
| Long solves blocking the page | frozen UI | Web Worker, `time_limit`, hard Worker timeout, progress text "homotopy …" |
| Users trying inputs outside v0.1 (intermediate/compound waves, Bt → 0 on both sides, switch-on/off, 180° problems) | confusion | show `reason` in plain language (table of messages per reason), link to the report |
| Numerical dependencies change under the reference project | golden data not reproducible | committed `reference/Manifest.toml`; golden files record commit and package versions; CI test compares them |
| Missing scripts for the small-Bt and 3,800-problem validations | golden data that don't match the README claims | write the scripts in Phase 1, check them against the README figures, ask the Julia author if they disagree |

---

## 8. Useful resources

- Julia source at the snapshot: [ExactMHDRiemannSolver](https://github.com/torrilhon/ExactMHDRiemannSolver), especially `src/`, `test/runtests.jl`, `test/data/` and `benchmark/stress.jl`.
- Vern9 implementation and tableaus: [OrdinaryDiffEqVerner sources](https://github.com/SciML/OrdinaryDiffEq.jl/tree/master/lib/OrdinaryDiffEqVerner/src) and [docs](https://sciml.github.io/OrdinaryDiffEq.jl/stable/explicit/Verner/).
- Verner's own coefficient files: [Jim Verner's Refuge for Runge–Kutta Pairs](https://www.sfu.ca/~jverner/), and the [paper on optimal pairs with interpolants](https://link.springer.com/article/10.1007/s11075-009-9290-3).
- DOP853 fallback: [Hairer's Fortran codes](https://www.unige.ch/~hairer/software.html), [modern Fortran edition](https://github.com/jacobwilliams/dop853).
- Internal numerical differentiation (frozen step sequences): H. G. Bock, "Numerical treatment of inverse problems in chemical reaction kinetics", Springer Series in Chemical Physics 18 (1981).
- Trust region: [NonlinearSolve.jl solver docs](https://docs.sciml.ai/NonlinearSolve/stable/native/solvers/) and [paper](https://arxiv.org/html/2403.16341v3); [MINPACK hybrj](https://www.math.utah.edu/software/minpack/minpack/hybrj.html) ([netlib](https://www.netlib.org/minpack/)) as a classical reference.
- Eigenvalues in JS: [ml-matrix](https://github.com/mljs/matrix) ([`EigenvalueDecomposition`](https://mljs.github.io/matrix/classes/EigenvalueDecomposition.html)).
- Brent: Roots.jl source (`Brent` in the bracketing methods) as the reference; JS cross-checks: [gist 1](https://gist.github.com/ryanspradlin/18c1010b7dd2d875284933d018c5c908), [gist 2](https://gist.github.com/borgar/3317728).
- Julia number printing and ranges: `Base.Ryu` (`repr` of Float64) and `Base.TwicePrecision` / `StepRangeLen` in the Julia sources.
- Testing: [Vitest `expect`](https://vitest.dev/api/expect.html). Use relative-error helpers rather than `toBeCloseTo`, which only counts absolute decimal digits.

---

## 9. Decisions

Settled:

1. Repository: [torrilhon/ExactMHDRiemann-Web](https://github.com/torrilhon/ExactMHDRiemann-Web).
2. License: MIT, like the original (already in the repository; the borrowed coefficient and algorithm sources are MIT as well).
3. Order: regular MHD solver first, with the report example as the first end-to-end target; quasi-Euler afterwards (§6).

Open:

1. Custom domain for the Pages URL later, and whether the repository stays public (Pages on a private repository needs a paid plan).
2. Plot library: uPlot (light) vs Plotly (built-in PNG export).
3. Vern9 vs DOP853 as the fan integrator: decided in Phase 2 by which one is verified first.
4. Whether the downloaded JSON should include the full dense fan data, so a later interactive version can re-sample without solving again. Probably not needed for v0.1.
