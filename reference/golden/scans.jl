# Sets 5 and 9: the validation scans described in the README of the Julia package, whose
# scripts are not in the Julia repository. These are reconstructions from the README
# text; their summaries are printed and stored so they can be compared with the README.

using DelimitedFiles, Random

# Set 5, README: "Small transverse field, 30 random problems per setting, both regimes
# c_A ≷ a on the small side, larger side Bt/√p from 1e-2 down to 1e-4 and smaller side
# from 1e-4 down to 1e-6: for a > c_A at least 29/30 Success in every setting; for
# c_A > a 27/30 (3 genuine vacuum limits) while the larger side is ≥ 1e-3, and 24-25/30
# when it is 3e-4 or 1e-4, the rest ending as NoConvergence."
#
# Reconstruction: problems 1-30 of random_problems.csv; the left state gets the larger
# Bt/√p, the right state the smaller one (field direction kept); the pressures of both
# states are scaled by one factor so that on the small side c_A = 2a (c_A > a) or
# a = 2c_A (a > c_A).
function gen_small_bt_scan()
    P = readdlm(joinpath(DATA, "random_problems.csv"), ',')
    γ = 5 / 3
    setbt(W, b) = (W = copy(W); f = b * sqrt(W[8]) / hypot(W[6], W[7]); W[6] *= f; W[7] *= f; W)
    larger = (1e-2, 3e-3, 1e-3, 3e-4, 1e-4)
    smaller = (1e-4, 1e-5, 1e-6)
    cases = Any[]
    summary = Any[]
    for regime in ("a>cA", "cA>a"), bl in larger, bs in smaller
        bs <= bl || continue
        sub = Any[]
        for i in 1:30
            L, R = copy(P[i, 2:9]), copy(P[i, 10:17])
            f = regime == "cA>a" ? (R[5]^2 / (4γ)) / R[8] : 4R[5]^2 / (γ * R[8])
            L[8] *= f; R[8] *= f
            L, R = setbt(L, bl), setbt(R, bs)
            rec = merge((regime = regime, larger = bl, smaller = bs, row = i, problem = probrec(L, R, γ)),
                first(golden_solve(L, R, γ; full = false)))
            push!(sub, rec); push!(cases, rec)
        end
        cnt = countmap_rc(sub)
        worst = maximum((c.check.maxerr for c in sub if c.retcode == Success); init = 0.0)
        push!(summary, (regime = regime, larger = bl, smaller = bs, counts = cnt, worst_maxerr = worst))
        println("  small-Bt $regime larger=$bl smaller=$bs: ", cnt, " worst maxerr ", worst)
    end
    save("solve_small_bt", (summary = summary, cases = cases))
end

# Set 9, README: "Bn = 0: Sod's problem matches Toro (p* = 0.30313, u* = 0.92745); 3800
# random perpendicular problems with Bt/√p from 0 to 1e3 and ratios up to 1:10⁴: all
# Success except cases that generate vacuum."
#
# Reconstruction: 19 levels of Bt/√p (0 and 18 log-spaced values in [1e-3, 1e3]) × 200
# problems. Per side ρ, p = exp(±log(1e4)/2 · U(-1,1)), so ratios up to 1:10⁴; vx in
# [-2, 2], vy, vz in [-0.5, 0.5]; Bt/√p = level · U(0.5, 1) in a random direction; Bn = 0;
# γ ∈ {1.4, 5/3, 2}. Plus Sod, the 200 problems of runtests.jl and the colliding flows.
function gen_solve_perp()
    cases = Any[]
    add!(name, L, R, γ) = push!(cases, merge((name = name, problem = probrec(L, R, γ)), first(golden_solve(L, R, γ))))
    add!("sod", [1.0, 0, 0, 0, 0, 0, 0, 1], [0.125, 0, 0, 0, 0, 0, 0, 0.1], 1.4)
    add!("sod_with_bt", [1.0, 0, 0, 0, 0, 1, 0, 1], [0.125, 0, 0, 0, 0, 0.5, 0, 0.1], 5 / 3)
    add!("colliding_1e8", [1.0, 1e8, 0, 0, 0, 1, 0, 1], [1.0, -1e8, 0, 0, 0, 1, 0, 1], 5 / 3)
    add!("diverging_vacuum", [1.0, -10, 0, 0, 0, 0.3, 0, 1], [1.0, 10, 0, 0, 0, 0.3, 0, 1], 5 / 3)
    rng = MersenneTwister(5)                         # runtests.jl, "quasi-Euler solver"
    for k in 1:200
        st() = [0.2 + 3rand(rng), 2rand(rng) - 1, rand(rng) - 0.5, rand(rng) - 0.5, 0.0,
                rand(rng, (0.0, 1e-8, 1.0)) * (2rand(rng) - 1), rand(rng, (0.0, 1.0)) * (2rand(rng) - 1), 0.2 + 3rand(rng)]
        L, R = st(), st()
        add!("runtests/$k", L, R, rand(rng, (1.4, 5 / 3, 2.0)))
    end
    levels = [0.0; 10.0 .^ range(-3, 3; length = 18)]
    rng = MersenneTwister(16)
    scan = Any[]
    summary = Any[]
    for b in levels
        sub = Any[]
        for _ in 1:200
            function side()
                ρ = exp((2rand(rng) - 1) * log(1e4) / 2); p = exp((2rand(rng) - 1) * log(1e4) / 2)
                bt = b * (0.5 + 0.5rand(rng)) * sqrt(p); th = 2π * rand(rng)
                [ρ, 4rand(rng) - 2, rand(rng) - 0.5, rand(rng) - 0.5, 0.0, bt * cos(th), bt * sin(th), p]
            end
            L, R = side(), side()
            γ = rand(rng, (1.4, 5 / 3, 2.0))
            rec = merge((level = b, problem = probrec(L, R, γ)), first(golden_solve(L, R, γ; full = false)))
            push!(sub, rec); push!(scan, rec)
        end
        cnt = countmap_rc(sub)
        push!(summary, (level = b, counts = cnt))
        println("  perpendicular Bt/√p = $b: ", cnt)
    end
    println("  perpendicular total: ", countmap_rc(scan))
    save("solve_perp", (cases = cases, scan_summary = summary, scan = scan))
end
