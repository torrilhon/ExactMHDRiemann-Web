# L2: the residual of the 5×5 system at many Ψ, including domain edges (where Julia
# returns the BIG vector) and the branch switches of build_side (ψf = 0, ψs = SLOW_EPS),
# plus ForwardDiff Jacobians at a subset of the points (plan §4.3, §5.1).

using ForwardDiff, DelimitedFiles

function canonical_problem(L, R, γ; opts = SolverOptions())
    F = E.make_frame(SVector{8}(L), SVector{8}(R))
    Lc, Rc = E.canonical_states(SVector{8}(L), SVector{8}(R), F)
    UL, UR = E.to_hstate(Lc), E.to_hstate(Rc)
    ctx = E.Ctx(γ, Lc[5], opts, E.problem_scale(UL, UR, Lc[5], γ))
    return UL, UR, ctx
end

function gen_residual()
    rng = MersenneTwister(15)
    P = readdlm(joinpath(DATA, "random_problems.csv"), ',')
    probs = Any[([3.0, 0, 0, 0, 1.5, 1.0, 0, 3.0], [1.0, 0, 0, 0, 1.5, cos(1.5), sin(1.5), 1.0], 5 / 3)]
    for i in 1:24
        push!(probs, (P[i, 2:9], P[i, 10:17], 5 / 3))
    end
    out = Any[]
    for (L, R, γ) in probs
        sol = solve(RiemannProblem(L, R; γ))
        UL, UR, ctx = canonical_problem(L, R, γ)
        Ψs = SVector{5,Float64}(sol.Ψ)
        pts = SVector{5,Float64}[Ψs]
        for sc in (1e-3, 0.1, 1.0), _ in 1:4
            push!(pts, Ψs + sc * SVector{5}(2rand(rng, 5) .- 1))
        end
        for _ in 1:8
            push!(pts, SVector(4rand(rng) - 2, 4rand(rng) - 2, 2π * rand(rng) - π, 4rand(rng) - 2, 4rand(rng) - 2))
        end
        # domain edges: saturated slow shocks, strong fans towards vacuum, huge shocks
        for e in (SVector(0.0, 50.0, 0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0, 200.0, 0.0),
                  SVector(-40.0, 0.0, 0.0, 0.0, 0.0), SVector(0.0, -60.0, 0.0, -60.0, 0.0), SVector(0.0, 0.0, 0.0, 0.0, 30.0))
            push!(pts, Ψs .* (e .== 0) + e)
        end
        # branch switches
        for k in (1, 5), v in (0.0, 1e-8, -1e-8, 1e-6, -1e-6)
            push!(pts, setindex(Ψs, v, k))
        end
        for k in (2, 4), v in (E.SLOW_EPS, E.SLOW_EPS + 1e-12, E.SLOW_EPS - 1e-12, 0.0, 1e-6, -1e-6)
            push!(pts, setindex(Ψs, v, k))
        end
        # the parametrizations keep the residual defined on most of R^5; domain errors
        # (DomainError, failed fan ODEs) appear only for |Ψ| of order 1e4: keep a few
        nbig = 0
        for _ in 1:400
            Ψ = SVector{5}(1e4 .* (2rand(rng, 5) .- 1))
            all(==(E.BIG), E.safe_residual(Ψ, UL, UR, ctx)) || continue
            push!(pts, Ψ); nbig += 1
            nbig == 4 && break
        end
        for _ in 1:4
            push!(pts, SVector{5}(1e4 .* (2rand(rng, 5) .- 1)))
        end
        evals = map(pts) do Ψ
            r = E.safe_residual(Ψ, UL, UR, ctx)
            (psi = collect(Ψ), r = collect(r), big = all(==(E.BIG), r))
        end
        # Jacobians: at the solution, two nearby points and one switch point
        jpts = [Ψs, pts[2], pts[6], setindex(Ψs, 1e-6, 1)]
        jacs = map(jpts) do Ψ
            f(x) = E.safe_residual(x, UL, UR, ctx)
            (psi = collect(Ψ), r = collect(f(Ψ)), J = Matrix(ForwardDiff.jacobian(f, Ψ)))
        end
        push!(out, (problem = probrec(L, R, γ), retcode = sol.retcode, UL = ser(UL), UR = ser(UR), ctx = ser(ctx),
            evals = evals, jacobians = jacs))
    end
    save("residual", (BIG = E.BIG, SLOW_EPS = E.SLOW_EPS, problems = out))
end
