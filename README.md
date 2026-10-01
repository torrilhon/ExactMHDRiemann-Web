# Exact MHD Riemann Solver for the web

Exact solutions of Riemann problems of one-dimensional ideal magnetohydrodynamics
(γ-law gas) in the browser. This is a TypeScript port of the Julia package
[ExactMHDRiemannSolver.jl](https://github.com/torrilhon/ExactMHDRiemannSolver), taken
as a snapshot of commit `6da1c02` (version 0.1.0), following

> M. Torrilhon, *Exact Solver and Uniqueness Conditions for Riemann Problems of Ideal
> Magnetohydrodynamics*, SAM Research Report 2002-06, ETH Zürich.

**Status:** under construction. See [PORT_PLAN.md](PORT_PLAN.md) for the plan and its
phases.

## Development

```sh
npm ci
npm run dev        # dev server
npm run check      # typecheck, numerics lint, tests, production build
```

The Julia project in [`reference/`](reference/) generates the golden test data; it is
not needed for building or testing.

## License

MIT, see [LICENSE](LICENSE).
