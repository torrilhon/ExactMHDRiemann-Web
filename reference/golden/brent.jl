# Brent (Roots.jl) on polynomials: roots and evaluation counts, for a bit-level test of
# the TS port of the method (plan §4.4). Polynomials only, so no libm differences.

using Roots

function gen_brent()
    rng = MersenneTwister(17)
    cases = Any[]
    function run(name, c, a, b; kw...)
        n = Ref(0)
        f(x) = (n[] += 1; evalpoly(x, Tuple(c)))
        res = try
            find_zero(f, (a, b), Roots.Brent(); kw...)
        catch err
            err isa Roots.ConvergenceFailed ? "ConvergenceFailed" :
            err isa ArgumentError ? "ArgumentError" : rethrow()
        end
        push!(cases, (name = name, c = collect(Float64, c), a = a, b = b, opts = Dict(string(k) => v for (k, v) in kw),
            root = res, evals = n[]))
    end
    # random cubics and quintics with a root in the bracket, the solver's tolerances
    for k in 1:150
        r = 4rand(rng) - 2
        deg = rand(rng, (1, 3, 5))
        c = randn(rng, deg + 1)
        c[1] -= evalpoly(r, Tuple(c))                     # root at r
        a, b = r - 3rand(rng), r + 3rand(rng)
        sign(evalpoly(a, Tuple(c))) * sign(evalpoly(b, Tuple(c))) < 0 || continue
        run("poly/$k/xrtol4eps", c, a, b; xatol = 0.0, xrtol = 4eps())
        run("poly/$k/sample", c, a, b; xatol = 1e-14)
        run("poly/$k/default", c, a, b)
    end
    # wide brackets as in perp_fan: root of ρ^γ-like power law replaced by polynomial in ρ
    for r in (1e-8, 1e-3, 0.3, 0.9999)
        run("wide/$r", [-r^2, 0.0, 1.0], 1e-300, 1.0; xatol = 0.0, xrtol = 4eps())
        run("wide_lin/$r", [-r, 1.0], 1e-300, 1.0; xatol = 0.0, xrtol = 4eps())
    end
    # exact zeros at an endpoint, no sign change, double root, steep and flat functions
    run("zero_at_b", [-1.0, 1.0], 0.0, 1.0; xatol = 0.0, xrtol = 4eps())
    run("zero_at_a", [-1.0, 1.0], 1.0, 2.0; xatol = 0.0, xrtol = 4eps())
    run("no_sign_change", [1.0, 0.0, 1.0], -1.0, 1.0; xatol = 0.0, xrtol = 4eps())
    run("reversed", [-0.5, 1.0], 1.0, 0.0; xatol = 0.0, xrtol = 4eps())
    run("flat_cubic", [0.0, 0.0, 0.0, 1.0], -1.0, 2.0; xatol = 0.0, xrtol = 4eps())
    run("steep", [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1e10], 0.0, 1.0; xatol = 0.0, xrtol = 4eps())
    run("tiny_root", [-1e-200, 1.0], -1.0, 1.0; xatol = 0.0, xrtol = 4eps())
    run("maxevals_hit", [-1e-200, 1.0], -1.0, 1.0; xatol = 0.0, xrtol = 0.0)
    save("brent", (cases = cases,))
end
