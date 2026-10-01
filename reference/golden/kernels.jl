# L1 kernels: shocks, rotations, slow limit, fan end states, the fast-shock cubic and the
# quasi-Euler kernels. Random states as in test/runtests.jl.

randstate(rng; Bn = 0.3 + 2rand(rng)) = E.HState(0.2 + 3rand(rng), 2rand(rng) - 1, 0.2 + 3rand(rng),
    0.1 + 2rand(rng), 2π * rand(rng), SVector(rand(rng) - 0.5, rand(rng) - 0.5)), Bn

# The preamble of E.fast_shock up to the root of the cubic (same formulas, Float64 only),
# so the TS port can test largest_root_in on exactly the coefficients the solver sees.
function cubic_case(U, ψ, ctx)
    γ, κ, Bn = ctx.γ, ctx.κ, ctx.Bn
    sp = sqrt(U.p)
    A, B = U.bt / sp, Bn / sp
    a0 = sqrt(γ * U.p / U.ρ)
    cf, cA, cs = E.speeds(U, ctx)
    df2 = cA > a0 ? cA^2 * U.bt^2 / (U.ρ * (cA^2 - cs^2)) : cf^2 - cA^2
    Df = γ * df2 / a0^2
    D = (sqrt(Df) + ψ)^2
    X = B^2 + D
    c = E.fast_cubic_bt(A, B^2, X, D, κ)
    yhi = sqrt(A^2 + 2 * (1 + X)) * (1 + 1e-12)
    y0 = E.largest_root_in(c, 1e-300, yhi)
    return (c = collect(c), lo = 1e-300, hi = yhi, root = y0)
end

function fast_rec(U, ψ, σ, ctx)
    D, s = E.fast_shock(U, ψ, σ, ctx)
    return (psi = ψ, sigma = σ, down = ser(D), speed = s, cubic = cubic_case(U, ψ, ctx))
end

ctxrec(ctx) = (gamma = ctx.γ, Bn = ctx.Bn)

function gen_kernels()
    # shocks, slow limit, rotation (runtests "kernels: Rankine-Hugoniot")
    rng = MersenneTwister(1)
    shocks = Any[]
    for γ in (5 / 3, 1.4, 2.0), _ in 1:40
        U, Bn = randstate(rng)
        ctx = E.Ctx(γ, Bn, SolverOptions())
        σ = rand(rng, (-1, 1))
        fast = [fast_rec(U, ψ, σ, ctx) for ψ in (1e-10, 1e-4, 0.05, 0.5, 3.0)]
        Δmax = E.slow_limit(U, ctx)
        slow = map((1e-7, 1e-3, 0.3, 0.9, 0.999)) do f
            D, s, M, v = E.slow_shock_state(U, f * Δmax, σ, ctx)
            (delta = f * Δmax, down = ser(D), speed = s, M = M, v = v,
                volume = collect(E.slow_volume(f * Δmax, U.bt / sqrt(U.p), (Bn / sqrt(U.p))^2, ctx.κ)),
                defect = E.slow_defect(U, f * Δmax, ctx))
        end
        Dr, sr = E.rotation(U, U.φ + 1.3, σ, ctx)
        push!(shocks, (ctx = ctxrec(ctx), U = ser(U), sigma = σ, fast = fast, slow_limit = Δmax, slow = slow,
            rotation = (alpha = U.φ + 1.3, down = ser(Dr), speed = sr)))
    end

    # fan end states and the checker's eigenvector integration (runtests "fans vs eigenvector")
    rng = MersenneTwister(2)
    fans = Any[]
    for γ in (5 / 3, 2.0), _ in 1:15
        U, Bn = randstate(rng)
        ctx = E.Ctx(γ, Bn, SolverOptions())
        W0 = E.to_prim(U, ctx)
        for σ in (-1, 1)
            Df, _ = E.fast_fan(U, -0.7, σ, ctx)
            Ds, _ = E.slow_fan(U, -0.4, σ, ctx)
            push!(fans, (ctx = ctxrec(ctx), U = ser(U), sigma = σ,
                fast = (psi = -0.7, down = ser(Df), eigen = collect(E.fan_eigen_integrate(W0, Df.ρ, σ < 0 ? 1 : 7, γ))),
                slow = (psi = -0.4, send = E.slow_fan_send(-0.4, ctx), down = ser(Ds),
                    eigen = collect(E.fan_eigen_integrate(W0, Ds.ρ, σ < 0 ? 3 : 5, γ)))))
        end
    end

    # near switch-on: tiny upstream Bt with cA > a (runtests "small transverse field")
    switch_on = Any[]
    ctx = E.Ctx(5 / 3, 2.0, SolverOptions())
    for A in (1e-3, 1e-5, 1e-7), ψ in (1e-8, 1e-4, 0.1, 1.0)
        U = E.HState(1.0, 0.0, 1.0, A, 0.0, SVector(0.0, 0.0))
        push!(switch_on, (ctx = ctxrec(ctx), U = ser(U), fast = fast_rec(U, ψ, 1, ctx)))
    end
    # slow limit on tiny Bt, both regimes (a > cA and cA > a)
    slow_small = Any[]
    for Bn in (0.5, 2.0), A in (1e-2, 1e-4, 1e-6)
        ctx = E.Ctx(5 / 3, Bn, SolverOptions())
        U = E.HState(1.0, 0.0, 1.0, A, 0.3, SVector(0.0, 0.0))
        Δmax = E.slow_limit(U, ctx)
        push!(slow_small, (ctx = ctxrec(ctx), U = ser(U), slow_limit = Δmax))
    end
    # vanishing normal field, regular kernels forced (runtests "small normal field")
    small_bn = Any[]
    for Bn in (1e-7, 1e-9, 1e-10)
        ctx = E.Ctx(5 / 3, Bn, SolverOptions(bn_euler = 0.0))
        U = E.HState(1.0, 0.0, 1.0, 0.8, 0.0, SVector(0.0, 0.0))
        for ψ in (1e-6, 0.1, 1.0)
            push!(small_bn, (ctx = ctxrec(ctx), U = ser(U), fast = fast_rec(U, ψ, 1, ctx)))
        end
    end

    # quasi-Euler kernels (Bn = 0): perp_side for P above and below the total pressure
    rng = MersenneTwister(13)
    perp = Any[]
    for _ in 1:60
        γ = rand(rng, (1.4, 5 / 3, 2.0))
        ctx = E.Ctx(γ, 0.0, SolverOptions())
        U = E.HState(0.2 + 3rand(rng), 2rand(rng) - 1, 0.2 + 3rand(rng), rand(rng, (0.0, 1e-8, 0.5, 2.0)) * rand(rng),
            2π * rand(rng), SVector(rand(rng) - 0.5, rand(rng) - 0.5))
        P0 = E.ptot(U.ρ, U.p, U.bt)
        σ = rand(rng, (-1, 1))
        for f in (1e-10, 1e-4, 0.01, 0.3, 0.9, 0.999, 1.001, 1.1, 3.0, 100.0, 1e6)
            D, w = E.perp_side(U, f * P0, σ, ctx; record = true)
            push!(perp, (ctx = ctxrec(ctx), U = ser(U), sigma = σ, P = f * P0, down = ser(D),
                wave = w === nothing ? nothing : (kind = w.kind, s_left = w.s_left, s_right = w.s_right,
                    send = w.fan === nothing ? nothing : w.fan.par)))
        end
    end
    # the strong-shock refusal: P/P0 beyond ~1e14 throws DomainError
    ctx = E.Ctx(5 / 3, 0.0, SolverOptions())
    U = E.HState(1.0, 0.0, 1.0, 0.5, 0.0, SVector(0.0, 0.0))
    perp_refusal = map((1e12, 1e13, 1e14, 1e15, 1e16)) do f
        P = f * E.ptot(U.ρ, U.p, U.bt)
        res = try
            D, _ = E.perp_shock(U, P, 1, ctx); ser(D)
        catch err
            err isa DomainError || rethrow(); "DomainError"
        end
        (ctx = ctxrec(ctx), U = ser(U), P = P, result = res)
    end

    save("kernels", (shocks = shocks, fans = fans, switch_on = switch_on, slow_small = slow_small,
        small_bn = small_bn, perp = perp, perp_refusal = perp_refusal))
end

# L1 dense: fan interiors from the dense ODE output.
function gen_fans_dense()
    rng = MersenneTwister(14)
    τs = [0.0; sort(rand(rng, 48)); 1.0]
    rec(f, ctx) = (family = f.family, sigma = f.σ, up = ser(f.up), par = f.par,
        tau = τs, states = [ser(E.fan_state(f, τ, ctx)) for τ in τs], speeds = [E.fan_speed(f, τ, ctx) for τ in τs])
    fans = Any[]
    for i in 1:100
        U, Bn = randstate(rng)
        γ = rand(rng, (1.4, 5 / 3, 2.0))
        ctx = E.Ctx(γ, Bn, SolverOptions())
        σ = rand(rng, (-1, 1))
        if isodd(i)
            ψ = -3rand(rng)
            _, f = E.fast_fan(U, ψ, σ, ctx; record = true)
        else
            ψ = -2rand(rng)
            _, f = E.slow_fan(U, ψ, σ, ctx; record = true)
        end
        push!(fans, merge((ctx = ctxrec(ctx), psi = ψ), rec(f, ctx)))
    end
    for _ in 1:20                                   # quasi-Euler fans (:fast0)
        γ = rand(rng, (1.4, 5 / 3, 2.0))
        ctx = E.Ctx(γ, 0.0, SolverOptions())
        U = E.HState(0.2 + 3rand(rng), 2rand(rng) - 1, 0.2 + 3rand(rng), 2rand(rng), 0.0, SVector(0.0, 0.0))
        P = E.ptot(U.ρ, U.p, U.bt) * (0.01 + 0.9rand(rng))
        σ = rand(rng, (-1, 1))
        _, f = E.perp_fan(U, P, σ, ctx; record = true)
        push!(fans, merge((ctx = ctxrec(ctx), P = P), rec(f, ctx)))
    end
    save("fans_dense", (fans = fans,))
end
