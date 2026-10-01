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

## Golden data generation (Phase 1)

All solve sets run with `time_limit = Inf`. Summary of the 8021 solves:

| set | problems | result | worst `check.maxerr` (Success) |
| --- | --- | --- | --- |
| examples and edge cases | 117 | 98 Success, 2 trivial, 15 refusals as intended, 2 RegularLimit, **2 CheckFailed** (below) | 2.1e-10 |
| random (`random_problems.csv`) | 2000 | 2000 Success (all direct) | 8.7e-9 |
| stress | 1000 | 978 Success, 22 RegularLimit vacuum | 1.95e-8 |
| small-Bt scan (reconstructed) | 900 | 835 Success, 45 RegularLimit vacuum, 20 NoConvergence | 1.1e-9 |
| quasi-Euler (incl. 3800 scan) | 4004 | 3925 Success, 78 RegularLimit vacuum, 1 NoConvergence (colliding flows 1e8, intended) | 6.0e-10 |

### CheckFailed on valid input (Julia bug, fixed on a branch)

`bt(1e-4) | bt(1e-6)` and `bt(2e-4) | bt(2e-6)` with `bt(b) = [1, 0, 0, 0, 1, b, 0, 1]`
are inside the supported domain (the second is in `runtests.jl`, which only tests
`retcode != Unsupported`) and return `CheckFailed`: "fast_fan: fan mismatch 6.1e-6"
(1.6e-6). The solution is right; its left fast fan is vanishingly weak (ψ ≈ −1e-8), and
`fan_pointwise_defect` divides by |dq/dτ| ~ 1e-8, so dense-output noise (≈ 8e-14/|ψ|)
is reported as a mismatch above `100·check_tol`.

Fix: commit `0ca567a` on branch `claude/elegant-edison-fb6a3w` of the Julia repository
(floor of 1e-6 of the state size in the denominator, plus a regression test; test suite
3289/3289). Regenerating all solve sets with the fix changes exactly these two cases to
Success; every other retcode, Ψ, wave table and check error is bit-identical. The golden
data stay at `6da1c02` until the fix is merged and the snapshot re-pinned.

### Time limit

Five homotopy solves of the small-Bt scan take 6.8–8.3 s in Julia, above the default
`time_limit = 5 s`; with the default they would end as NoConvergence. The
random problems include direct solves of up to 1.2 s (problem 42). The TS port will be
slower, so the shipped `time_limit` must come from measured TS times (plan §4.2).

### Domain errors in the residual

The parametrizations (tanh caps of slow shocks and slow fans) keep the residual defined
on most of R^5: no `BIG` vector in boxes up to |Ψ| ≤ 50 on any tested problem. Domain
errors (DomainError, fan ODEs failing with NaN step sizes) appear for |Ψ| of order 1e4.

### Reconstructed scans versus the README

- Small-Bt, c_A > a: 27/30 (3 vacuum) for the larger side ≥ 1e-3, 25/30 at 3e-4, 24/30
  at 1e-4: **as in the README**. a > c_A: 30/30, except 28–29/30 with the larger side at
  1e-4 (README: at least 29/30). The regime is set by scaling both pressures so that
  a = 2c_A on the small side; the README's construction is unknown.
- Perpendicular: 3723 Success, 77 vacuum of 3800: **as in the README** ("all Success
  except cases that generate vacuum").

### Hardware dependence of the reference: FMA

Julia's `evalpoly` (used for the fast-shock cubic) is Horner with `muladd`, which Julia
compiles to a fused multiply-add where the CPU has one (x86-64 with FMA3, ARM64), and to
a separate multiply and add otherwise. The golden data were generated on a machine with
FMA, so they are reproduced bit for bit only with fma semantics. The TS port emulates
fma exactly (`src/solver/math.ts`, tested on 3012 triples against Julia). OrdinaryDiffEq's
Verner steps also use `@muladd`; the ODE results are compared with tolerances anyway.

### FMA inside the solver, libm, and coplanar problems (Phase 3)

The fma question is subtler than above. `muladd` only *permits* fusion; LLVM decides
per compilation context. Inside the solver's inlined, type-stable code it mostly does
not fuse: plain multiply-add reproduces Julia's fast-shock cubic roots bit for bit in
600 of 621 golden cases, fused fma in 476 (the golden generators that call `evalpoly`
through dynamic dispatch get the fused version). The port therefore uses plain
arithmetic, and the remaining cases differ in the last bits (Brent's tolerance is 4 ULP).
The exact fma emulation is kept in the tests (`test/fma.ts`), where it reproduces the
Brent and fma golden data bit for bit.

`Math.sin`/`Math.cos` and Julia's differ in the last bit for some arguments. This
matters for **exactly coplanar problems in a rotated frame**: the canonical Bz of the
right state is mathematically 0, numerically ±1e-17, in Julia as in the port, so the
twist angle α = ±π and with it the first start (±π/2) is decided by round-off. The two
outcomes are mirror images (z → -z in the canonical frame) of the same regular solution
family. Axis-aligned coplanar problems (Brio–Wu: Bz = 0 exactly) are not affected. Golden
comparisons of coplanar problems accept the mirrored solution.

## Port vs reference, final comparison (Phases 5–6)

`npm run compare-golden <set>`, all with `time_limit = Inf`:

| set | identical retcode/reason/method | wave tables | TS time (mean / max) | Julia total |
| --- | --- | --- | --- | --- |
| examples and edge cases | 115/117 (the 2 fixed CheckFailed) | ≤ 4.1e-13 | 13 ms / 0.7 s | 0.8 s |
| random | 2000/2000 | ≤ 1.8e-14 | 3.4 ms / 36 ms | 4.1 s |
| stress | 1000/1000 | ≤ 2.2e-12 | 22 ms / 1.3 s | 9.7 s |
| small-Bt scan | 894/900 (99.3 %) | ≤ 1.4e-12 | 231 ms / 28 s | 67.8 s |
| quasi-Euler | 4004/4004 | ≤ 7.8e-13 | 1.5 ms / 21 ms | 3.8 s |

What made the difference (Phase 6): for static arrays, NonlinearSolve's
`MoreTrustRegionDescent` solves the **normal equations** with StaticArrays' LU, which
throws `SingularException` on an exactly zero pivot; `nl_solve` then gives up on that
start at once. In the small-Bt problems the Jacobian at the start vectors has an exactly
zero column (the saturated slow-shock variable of the small-Bt side), so Julia moves on
to the homotopy after microseconds. A port that regularizes the rank deficiency instead
(MINPACK `lmpar`, the first version) iterates on and visits very expensive fans (34 653
Vern9 steps each, in Julia too): up to 260 s per problem instead of 50 ms. The port now
follows the normal-equation path, including the exception.

The 6 remaining small-Bt differences (a > c_A, larger side Bt/√p ≤ 3e-4): 5 problems
that Julia solves by continuation end as NoConvergence in the port, and 1 the other way.
In these the homotopy starts at the degenerate point R = L, where the twist direction is
singular; Julia's successful paths let the twist angle wind through ~30 turns. Finite
differences and ForwardDiff cannot be expected to agree there. Success counts: Julia 835,
port 831 of 900.

## Re-pin to `0ca567a`

The checker fix was merged into the Julia `main` (`0ca567ab6446e7da404a4016df8d0f484ae30595`)
and the reference project re-pinned to it; only the ExactMHDRiemannSolver entry of the
Manifest changed. Julia test suite: 3289/3289. All golden data were regenerated; the
two former CheckFailed cases are now Success in the reference, and the port matches
the examples set 117/117. Sections above that name `6da1c02` describe the state before
the re-pin.

## Re-pin to the release v0.1.0 (`0e8fdb7`)

The owner published the Julia release
[v0.1.0](https://github.com/torrilhon/ExactMHDRiemannSolver/releases/tag/v0.1.0) at
`0e8fdb7690c47be797c999d7ef551d6e34586e30`. It differs from `0ca567a` only in
`README.md` and `CITATION.cff` (`src/` is identical), so the reference project was
re-pinned to it: only `repo-rev` and `git-tree-sha1` of the ExactMHDRiemannSolver entry
of the Manifest changed. All golden data were regenerated. Apart from the `commit` in
each `meta` entry and the measured `time_ms` of the solve sets, every file is identical
to the data at `0ca567a`.
