# Exact MHD Riemann Solver for the web

Exact solutions of Riemann problems of one-dimensional ideal magnetohydrodynamics
(γ-law gas), computed in the browser: **https://torrilhon.github.io/ExactMHDRiemann-Web/**

This is a TypeScript port of the Julia package
[ExactMHDRiemannSolver.jl](https://github.com/torrilhon/ExactMHDRiemannSolver), taken as
a snapshot of commit `6da1c02` (version 0.1.0), following

> M. Torrilhon, *Exact Solver and Uniqueness Conditions for Riemann Problems of Ideal
> Magnetohydrodynamics*, SAM Research Report 2002-06, ETH Zürich.

The page solves regular problems (fast and slow shocks and fans, rotational
discontinuities, contact) and, for a vanishing normal field, the quasi-Euler problem
(fast waves and a tangential discontinuity). Every solution is checked independently
(Rankine–Hugoniot, Lax and entropy conditions, eigenvector integration of the fans);
inputs outside the supported domain are refused with a reason. Results can be
downloaded as `solution.csv` (the format of Julia's `write_csv`), `solution.json`, and
`problem.toml`, which Julia's `problem_load` reads, so every result can be reproduced
with the Julia package.

## Accuracy

The port is tested against reference data generated once by the Julia package
([`test/golden/`](test/golden/README.md), [`reference/`](reference/README.md)):

| test set | problems | agreement with Julia |
| --- | --- | --- |
| report example, Brio–Wu, C-code benchmarks, refusals, thresholds, symmetries, small fields | 117 | same retcode in all cases but two (Julia's checker false positive, fixed); wave tables ≤ 4e-13 |
| random problems of the Julia test data | 2000 | all identical; ≤ 8e-15; C reference solutions ≤ 1.3e-8 (8 decimals) |
| stress set (`benchmark/stress.jl`) | 1000 | all identical; ≤ 2.2e-12 |
| small transverse field scan | 900 | 894 identical, 6 degenerate continuation cases differ (Julia 835 Success, port 831); ≤ 1.4e-12 |
| quasi-Euler (Sod, 3800-problem scan) | 4004 | all identical; ≤ 8e-13; Sod within 1e-5 of Toro |
| kernels, fans, residual, profiles, CSV | — | shocks ≤ 1e-15, fan interiors ≤ 8e-13, sampled profiles ≤ 2.4e-13 |

Known differences and how the comparison treats round-off ties are listed in
[`test/golden/DIFFERENCES.md`](test/golden/DIFFERENCES.md); the porting plan and its
findings are in [`PORT_PLAN.md`](PORT_PLAN.md).

## Development

```sh
npm ci
npm run dev                         # dev server
npm run check                       # typecheck, numerics lint, tests, production build
npm run fuzz -- 2000                # random problems, fails on any failed check
npm run compare-golden solve_stress # full comparison of a golden set with Julia
```

The Julia project in [`reference/`](reference/) generates the golden data; it is not
needed for building or testing.

## Citation and license

Please cite the report above ([`CITATION.cff`](CITATION.cff)). MIT license, see
[LICENSE](LICENSE).
