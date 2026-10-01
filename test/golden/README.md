# Golden data

Generated **once** by the Julia reference project (`reference/generate_golden.jl`, plan
§5.1) at ExactMHDRiemannSolver `0e8fdb7` (release v0.1.0) with the package versions of
`reference/Manifest.toml`; every file records both in its `meta` entry, and
`test/unit/golden-meta.test.ts` checks them. Do not edit by hand; regenerate with

```sh
julia --project=reference reference/generate_golden.jl [part ...]
```

Files are gzipped JSON (`test/golden.ts` loads them). Numbers are Julia `repr`
(round-trip exact); `Inf`, `-Inf`, `NaN` are stored as strings. All solves ran with
`time_limit = Inf` (plan §4.2); `time_ms` is Julia's wall time per solve, for the
performance comparison only (it changes on every regeneration).

| file | level | content |
| --- | --- | --- |
| `basics` | L0 | `speeds`, `conserved`, `flux`, `primitive_jacobian` for 500 states; canonical frame of 300 L/R pairs (frame, canonical states, HStates, round trip, `speed_to_user`) |
| `format` | L0 | Julia `repr` of ~550 Float64 values; `range(a, b; length = n)` for 14 grids |
| `kernels` | L1 | 120 states (runtests `MersenneTwister(1)`): `fast_shock` at 5 strengths with the cubic's coefficients, bracket and root; `slow_limit`; `slow_shock_state` at 5 fractions with `slow_volume` and `slow_defect`; `rotation`. 60 fan cases (`MersenneTwister(2)`): `fast_fan(-0.7)`, `slow_fan(-0.4)` end states and the checker's `fan_eigen_integrate`. Switch-on, small-Bt slow limit, small-Bn (`bn_euler = 0`) shocks. Quasi-Euler `perp_side` at 11 pressure ratios on 60 states, and the strong-shock refusal |
| `brent` | L1 | Roots.jl Brent on 455 polynomial cases: roots and evaluation counts |
| `trustregion` | L2 | NonlinearSolve `TrustRegion` on 16 classic test problems: retcode, solution, steps |
| `fans_dense` | L1 | 120 recorded fans (50 fast, 50 slow, 20 quasi-Euler `:fast0`): `fan_state` and `fan_speed` at 50 τ |
| `residual` | L2 | 25 problems × 56 Ψ: solution, perturbations, random box, branch switches (ψf = 0, ψs = `SLOW_EPS`), `|Ψ| ~ 1e4` (130 points give the `BIG` vector); ForwardDiff Jacobians at 4 points per problem |
| `solve_examples` | L3 | 117 problems: examples, C-code benchmarks, one case per refusal reason and the thresholds, trivial, symmetries, switch continuity, small Bn forced, vanishing waves, small Bt (sets 1, 6–8, 10–12); full records incl. canonical waves |
| `solve_random` | L3 | the 2000 problems of `data/random_problems.csv` (sets 2, 3) |
| `solve_stress` | L3 | the 1000 problems of `benchmark/stress.jl`, inputs exported (set 4) |
| `solve_small_bt` | L3 | small-Bt scan, 900 problems (set 5, reconstructed, see `reference/BASELINE.md`) |
| `solve_perp` | L3 | quasi-Euler: Sod, colliding/diverging flows, the 200 runtests problems, the 3800-problem scan (set 9, reconstructed) |
| `output` | L4 | `sample` on 401 ξ for 50 random problems; index of the CSV files |
| `csv/*.csv` | L4 | `write_csv` output of `paper.toml`, `briowu.toml` (their `[output]` grids) and Sod, for the character-for-character test |
| `data/` | — | `random_problems.csv`, `c_reference.csv` of the Julia package (test data, MIT) |
