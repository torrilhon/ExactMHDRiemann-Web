# Baseline of the reference project

Results of the Julia snapshot with the pinned `Manifest.toml`, recorded in Phase 0
(2026-10-01) to confirm that the pinned numerics reproduce the published figures.

| | |
| --- | --- |
| ExactMHDRiemannSolver | 0.1.0 @ `6da1c024b8e5af6975dd66569e38b936424cb0e9` |
| Julia | 1.11.7 |
| NonlinearSolve | 4.32.0 (NonlinearSolveFirstOrder 2.10.0, NonlinearSolveBase 2.54.1) |
| OrdinaryDiffEqVerner | 2.4.2 (OrdinaryDiffEqCore 4.18.1) |
| Roots | 3.0.9 |
| ForwardDiff | 1.4.6 |
| SciMLBase | 3.57.0 |
| StaticArrays | 1.9.22 |

## Test suite

`julia --project=reference test/runtests.jl`, run in a checkout of the Julia repository
at the snapshot commit: **3287 / 3287 pass** (31.7 s).

## Examples (`reference/smoke.jl`)

| problem | retcode | method | residual | check.maxerr |
| --- | --- | --- | --- | --- |
| `paper.toml` | Success | direct | 9.6e-17 | 1.8e-11 |
| `briowu.toml` | Success | direct | 5.6e-17 | 3.8e-11 |

The wave table of `paper.toml` agrees with Tables 1–2 of the report to the printed
digits (vz of the slow-fan tail 0.438329, as in the C code; the report prints 0.438321).

## Stress set (`benchmark/stress.jl 1000`)

| | |
| --- | --- |
| Success | 978 (97.8 %) |
| RegularLimit `:vacuum` | 22 (2.2 %) |
| NoConvergence / CheckFailed | 0 / 0 |
| worst `check.maxerr` among Success | **1.95e-8** |
| max time per solve | 0.37 s |

This matches the README. Note that the worst check error of a `Success` exceeds
`check_tol = 1e-8`: `check_waves` accepts fan mismatches up to `100·check_tol`, so
`maxerr` alone is not a pass/fail criterion (plan §5.3).
