# Smoke test of the reference project: solves the report example and Brio-Wu at the
# pinned snapshot and prints the wave tables. Run from the repository root:
#   julia --project=reference reference/smoke.jl
using ExactMHDRiemannSolver
using Pkg

const E = ExactMHDRiemannSolver
info = Pkg.dependencies()[Base.UUID("0d315d94-8d33-4e3b-9f90-006d7b2d27bf")]
println("ExactMHDRiemannSolver ", info.version, " @ ", info.git_revision, ", Julia ", VERSION)

for file in ("paper.toml", "briowu.toml")
    L, R, γ, t, x = problem_load(joinpath(pkgdir(E), "examples", file))
    t0 = time()
    sol = solve(RiemannProblem(L, R; γ))
    dt = time() - t0
    println("\n== ", file, ": ", sol.retcode, " method = ", sol.method, ", residual = ", sol.residual,
        ", check.maxerr = ", sol.check.maxerr, ", ", round(1e3dt, digits = 1), " ms (incl. compilation on first call)")
    show(stdout, MIME"text/plain"(), sol)
    println()
end
