# NonlinearSolve.TrustRegion (exactly as nl_solve calls it) on classic test problems
# (Moré, Garbow, Hillstrom, ACM TOMS 7, 1981), for the TS port of the trust region.

function tr_problems()
    P = Any[]
    push!(P, ("rosenbrock", u -> SVector(10 * (u[2] - u[1]^2), 1 - u[1]), [[-1.2, 1.0], [5.0, -3.0]]))
    push!(P, ("powell_badly_scaled", u -> SVector(1e4 * u[1] * u[2] - 1, exp(-u[1]) + exp(-u[2]) - 1.0001), [[0.0, 1.0]]))
    push!(P, ("helical_valley", u -> begin
        θ = atan(u[2], u[1]) / (2π) + (u[1] < 0 ? 0.5 : 0.0)
        SVector(10 * (u[3] - 10θ), 10 * (sqrt(u[1]^2 + u[2]^2) - 1), u[3])
    end, [[-1.0, 0.0, 0.0], [1.5, 1.0, 0.3]]))
    push!(P, ("powell_singular", u -> SVector(u[1] + 10u[2], sqrt(5) * (u[3] - u[4]), (u[2] - 2u[3])^2, sqrt(10) * (u[1] - u[4])^2),
        [[3.0, -1.0, 0.0, 1.0]]))
    push!(P, ("trigonometric5", u -> begin
        n = 5; s = sum(cos, u)
        SVector{5}([n - s + i * (1 - cos(u[i])) - sin(u[i]) for i in 1:n])
    end, [fill(0.2, 5), fill(1.0, 5)]))
    push!(P, ("broyden_tridiagonal5", u -> begin
        n = 5; x(i) = (i < 1 || i > n) ? 0.0 : u[i]
        SVector{5}([(3 - 2u[i]) * u[i] - x(i - 1) - 2x(i + 1) + 1 for i in 1:n])
    end, [fill(-1.0, 5), fill(3.0, 5)]))
    push!(P, ("singular_start", u -> SVector(u[1]^2 + u[2]^2 - 1, u[1] - u[2]), [[0.0, 0.0]]))
    push!(P, ("freudenstein_roth", u -> SVector(-13 + u[1] + ((5 - u[2]) * u[2] - 2) * u[2], -29 + u[1] + ((u[2] + 1) * u[2] - 14) * u[2]),
        [[0.5, -2.0], [6.0, 3.0]]))
    push!(P, ("linear5", u -> SVector{5}([sum((i == j ? 4.0 : 1.0 / (i + j)) * u[j] for j in 1:5) - i for i in 1:5]), [zeros(5)]))
    push!(P, ("big_plateau", u -> abs(u[1]) > 2 ? SVector(1e6, 1e6) : SVector(u[1]^3 - 0.5, u[2] - u[1]), [[1.9, 0.0], [3.0, 0.0]]))
    return P
end

function gen_trustregion()
    NS = E.NonlinearSolve
    cases = Any[]
    for (name, F, starts) in tr_problems(), u0 in starts
        n = length(u0)
        f(u, p) = F(u)
        prob = NS.NonlinearProblem{false}(f, SVector{n,Float64}(u0))
        s = NS.solve(prob, NS.TrustRegion(; autodiff = NS.AutoForwardDiff()); abstol = 1e-12, maxiters = 60)
        push!(cases, (name = name, u0 = u0, u = collect(s.u), retcode = string(s.retcode),
            resnorm = maximum(abs, F(s.u)), nsteps = s.stats.nsteps, nf = s.stats.nf, njacs = s.stats.njacs))
        println("  ", rpad(name, 22), u0, " → ", s.retcode, " ‖F‖∞ = ", maximum(abs, F(s.u)), " steps ", s.stats.nsteps)
    end
    save("trustregion", (cases = cases,))
end
