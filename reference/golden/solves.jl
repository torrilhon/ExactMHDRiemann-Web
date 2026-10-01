# L3: full solves (plan §5.2). All with time_limit = Inf (plan §4.2).

using DelimitedFiles, Random

const s4 = sqrt(4π)

caserec(name, L, R, γ; opts = SolverOptions(), full = true) =
    merge((name = name, problem = probrec(L, R, γ), opts = optsrec(opts)), first(golden_solve(L, R, γ; opts, full)))

# Sets 1, 6, 7, 8, 10, 11, 12 and the C-code benchmarks: small, full records.
function gen_solve_examples()
    cases = Any[]
    add!(name, L, R, γ = 5 / 3; kw...) = push!(cases, caserec(name, L, R, γ; kw...))
    P = readdlm(joinpath(DATA, "random_problems.csv"), ',')

    # set 1: examples
    for file in ("paper.toml", "briowu.toml")
        p = problem_load(joinpath(pkgdir(E), "examples", file))
        add!("examples/" * file, p.L, p.R, p.γ)
    end
    # C-code benchmarks of runtests.jl (Brio-Wu and two Ryu-Jones type problems)
    add!("c_benchmark/1", [1, 0, 0, 0, 0.75, 1, 0, 1], [0.125, 0, 0, 0, 0.75, -1, 0, 0.1], 2.0)
    add!("c_benchmark/2", [1.08, 1.2, 0.01, 0.5, 2 / s4, 3.6 / s4, 2 / s4, 0.95], [1, 0, 0, 0, 2 / s4, 4 / s4, 2 / s4, 1])
    add!("c_benchmark/3", [1, 10, 0, 0, 5 / s4, 5 / s4, 0, 20], [1, -10, 0, 0, 5 / s4, 5 / s4, 0, 1])

    # set 6: refusals, one per reason, and inputs at the thresholds
    good = [1.0, 0, 0, 0, 1.0, 1.0, 0, 1.0]
    bt(b) = [1.0, 0, 0, 0, 1.0, b, 0, 1.0]
    add!("refusal/nonfinite", [NaN; good[2:8]], good)
    add!("refusal/nonfinite_inf", good, [good[1:7]; Inf])
    add!("refusal/nonfinite_gamma", good, good, NaN)
    add!("refusal/gamma", good, good, 1.0)
    add!("refusal/density", [-1.0; good[2:8]], good)
    add!("refusal/density_zero", good, [0.0; good[2:8]])
    add!("refusal/pressure", good, [good[1:7]; 0.0])
    add!("refusal/bn_jump", good, [1.0, 0, 0, 0, 2.0, 1.0, 0, 1.0])
    add!("refusal/bn_jump_tol_ok", good, [1.0, 0, 0, 0, 1.0 + 0.9e-12, 1.0, 0, 1.0])
    add!("refusal/switch_on_off_rj4d", [1.0, 0, 0, 0, 0.7, 0, 0, 1.0], [0.3, 0, 0, 1, 0.7, 1, 0, 0.2])
    add!("refusal/extreme_ratio", good, [1e-7, 0, 0, 0, 1.0, 1.0, 0, 1.0])
    add!("refusal/extreme_ratio_bn0", [1.0, 0, 0, 0, 0, 1, 0, 1.0], [1e-7, 0, 0, 0, 0, 1, 0, 1.0])
    add!("refusal/extreme_ratio_at_limit", good, [1e-6, 0, 0, 0, 1.0, 1.0, 0, 1.0])
    add!("refusal/bt_smaller", good, bt(5e-7))
    add!("refusal/bt_both", bt(5e-5), bt(5e-5))
    add!("refusal/bt_smaller_at_threshold", good, bt(1e-6))
    add!("refusal/bt_larger_at_threshold", bt(1e-4), bt(1e-6))
    add!("refusal/bt_ok_1", bt(2e-6), good)
    add!("refusal/bt_ok_2", bt(2e-4), bt(2e-6))
    add!("refusal/bn_euler_at_threshold", [1.0, 0, 0, 0, 1e-10, 1, 0, 1.0], [0.125, 0, 0, 0, 1e-10, 0.5, 0, 0.1])
    add!("refusal/bn_euler_above", [1.0, 0, 0, 0, 1.0000001e-10, 1, 0, 1.0], [0.125, 0, 0, 0, 1.0000001e-10, 0.5, 0, 0.1])
    add!("refusal/vacuum_diverging", [1.0, -8, 0, 0, 1.0, 1.0, 0, 1.0], [1.0, 8, 0, 0, 1.0, 1.0, 0, 1.0])

    # set 7: trivial
    add!("trivial", good, good)
    add!("trivial_rotated", [1.0, 0.2, 0.1, 0, 1.0, 0.6, 0.8, 1.0], [1.0, 0.2, 0.1, 0, 1.0, 0.6, 0.8, 1.0])

    # set 8: symmetries (as in runtests.jl, MersenneTwister(3))
    rng = MersenneTwister(3)
    rot(W, t) = (c = cos(t); s = sin(t);
        [W[1], W[2], c * W[3] - s * W[4], s * W[3] + c * W[4], W[5], c * W[6] - s * W[7], s * W[6] + c * W[7], W[8]])
    negB(W) = [W[1:4]; -W[5:7]; W[8]]
    shift(W, u, v, w) = [W[1], W[2] + u, W[3] + v, W[4] + w, W[5:8]...]
    mirror(W) = [W[1], -W[2], W[3], W[4], -W[5], W[6], W[7], W[8]]
    for i in 1:8
        L = P[i, 2:9]; R = P[i, 10:17]; t = 2π * rand(rng)
        add!("symmetry/$i/base", L, R)
        add!("symmetry/$i/rotate", rot(L, t), rot(R, t))
        add!("symmetry/$i/negB", negB(L), negB(R))
        add!("symmetry/$i/shift", shift(L, 0.7, -0.3, 0.2), shift(R, 0.7, -0.3, 0.2))
        add!("symmetry/$i/mirror", mirror(R), mirror(L))
    end

    # set 10: switch continuity at bn_euler
    for i in 1:5
        L, R = copy(P[i, 2:9]), copy(P[i, 10:17]); b = sqrt(min(L[8], R[8]))
        L[5] = R[5] = 2e-10 * b; add!("switch/$i/regular", copy(L), copy(R))
        L[5] = R[5] = 0.9e-10 * b; add!("switch/$i/quasi_euler", copy(L), copy(R))
    end

    # set 11: small normal field with the regular solver forced
    o = SolverOptions(bn_euler = 0.0)
    for i in (15, 18, 27), bn in (1e-7, 1e-9, 1e-10)
        L, R = copy(P[i, 2:9]), copy(P[i, 10:17]); b = bn * sqrt(min(L[8], R[8])); L[5] = b; R[5] = b
        add!("small_bn/$i/$bn", L, R; opts = o)
    end

    # set 12: vanishing waves (single-wave problems)
    R = [1.0, 0, 0, 0, 0.8, 0.6, 0, 1.0]
    ctx = E.Ctx(5 / 3, 0.8, SolverOptions())
    U = E.to_hstate(R)
    for (kind, make) in ((:fast_shock, u -> E.fast_shock(u, 0.3, 1, ctx)[1]),
                         (:slow_shock, u -> E.slow_shock_state(u, 0.2, 1, ctx)[1]),
                         (:fast_fan, u -> E.fast_fan(u, -0.4, 1, ctx)[1]),
                         (:slow_fan, u -> E.slow_fan(u, -0.3, 1, ctx)[1]),
                         (:rotation, u -> E.rotation(u, 0.9, 1, ctx)[1]))
        L = Vector(E.to_prim(make(U), ctx))
        add!("single_wave/$kind", L, R)
    end

    # the small-Bt cases of runtests.jl (one side with Bt/√p = 1e-5, both regimes, both sides)
    setbt(W, b) = (W = copy(W); f = b * sqrt(W[8]) / hypot(W[6], W[7]); W[6] *= f; W[7] *= f; W)
    for i in 1:6, lowβ in (false, true)
        L, R = copy(P[i, 2:9]), copy(P[i, 10:17])
        if lowβ
            f = (R[5]^2 / (4 * 5 / 3)) / R[8]; L[8] *= f; R[8] *= f
        end
        add!("small_bt/$i/$(lowβ ? "cA>a" : "a>cA")/right", L, setbt(R, 1e-5))
        add!("small_bt/$i/$(lowβ ? "cA>a" : "a>cA")/left", setbt(R, 1e-5), L)
    end

    save("solve_examples", (cases = cases,))
end

# Sets 2 and 3: the 2000 random problems (compact records; c_reference.csv is copied
# to test/golden/data and compared in TS).
function gen_solve_random()
    P = readdlm(joinpath(DATA, "random_problems.csv"), ',')
    cases = map(1:size(P, 1)) do i
        L, R = P[i, 2:9], P[i, 10:17]
        merge((seed = Int(P[i, 1]), problem = probrec(L, R, 5 / 3)), first(golden_solve(L, R, 5 / 3; full = false)))
    end
    save("solve_random", (cases = cases,))
end

# Set 4: the stress set of benchmark/stress.jl, inputs exported.
function gen_solve_stress(n = 1000)
    function side(rng, Bn; bmax = 3.0, ratio = 100.0)
        bt = exp(rand(rng) * log(bmax / 0.01)) * 0.01; th = 2π * rand(rng)
        ρ = exp((2rand(rng) - 1) * log(ratio) / 2); p = exp((2rand(rng) - 1) * log(ratio) / 2)
        [ρ, 3 * (2rand(rng) - 1), 2rand(rng) - 1, 2rand(rng) - 1, Bn, bt * cos(th), bt * sin(th), p]
    end
    function gen(rng, i)
        Bn = exp((2rand(rng) - 1) * log(30)) * 0.3
        L = side(rng, Bn); R = side(rng, Bn)
        if i % 5 == 0
            R[6:7] = -hypot(R[6], R[7]) .* L[6:7] ./ hypot(L[6], L[7])
        end
        return L, R, rand(rng, (5 / 3, 1.4, 2.0))
    end
    rng = MersenneTwister(7)
    cases = map(1:n) do i
        L, R, γ = gen(rng, i)
        merge((index = i, problem = probrec(L, R, γ)), first(golden_solve(L, R, γ; full = false)))
    end
    summary = countmap_rc(cases)
    println("  stress: ", summary)
    save("solve_stress", (summary = summary, cases = cases))
end

function countmap_rc(cases)
    d = Dict{String,Int}()
    for c in cases
        k = string(c.retcode) * (c.reason === :none ? "" : ":" * string(c.reason))
        d[k] = get(d, k, 0) + 1
    end
    return d
end
