# Julia reference project

This Julia project produces the golden data in `test/golden/` that the TypeScript port
is tested against (plan §5.1). It is not needed to build or test the web app.

`Project.toml` depends on
[ExactMHDRiemannSolver](https://github.com/torrilhon/ExactMHDRiemannSolver) at the
snapshot commit `6da1c024b8e5af6975dd66569e38b936424cb0e9`. The committed
`Manifest.toml` pins **all** numerical dependencies (NonlinearSolve,
OrdinaryDiffEqVerner, Roots, ForwardDiff, …). The Julia repository has no Manifest of
its own, so this file, not the commit alone, fixes the reference numerics. A unit test
(`test/unit/snapshot.test.ts`) checks that the Manifest and `src/solver/snapshot.ts`
name the same commit.

```sh
julia --project=reference -e 'using Pkg; Pkg.instantiate()'
julia --project=reference reference/smoke.jl     # report example and Brio-Wu
```

Do not run `Pkg.update()` here: changing the Manifest changes the reference. If it is
ever needed, regenerate all golden data in the same commit.
