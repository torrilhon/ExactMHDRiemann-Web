# L4: sampled profiles on the write_csv grid, and the CSV files themselves for the
# character-for-character test of the TS writer (plan §4.7).

using DelimitedFiles

function gen_output()
    mkpath(joinpath(OUT, "csv"))
    samples = Any[]
    for file in ("paper.toml", "briowu.toml")
        p = problem_load(joinpath(pkgdir(E), "examples", file))
        sol = solve(RiemannProblem(p.L, p.R; γ = p.γ))
        name = splitext(file)[1]
        write_csv(joinpath(OUT, "csv", name * ".csv"), sol; t = p.t, x = p.x)
        push!(samples, (name = name, problem = probrec(p.L, p.R, p.γ), t = p.t, x = [first(p.x), last(p.x)],
            n = length(p.x), csv = "csv/$name.csv"))
    end
    # Sod (quasi-Euler) on the default grid
    sod = solve(RiemannProblem([1.0, 0, 0, 0, 0, 0, 0, 1], [0.125, 0, 0, 0, 0, 0, 0, 0.1]; γ = 1.4))
    write_csv(joinpath(OUT, "csv", "sod.csv"), sod; t = 0.2, x = range(-0.5, 0.5; length = 1001))
    push!(samples, (name = "sod", problem = probrec([1.0, 0, 0, 0, 0, 0, 0, 1], [0.125, 0, 0, 0, 0, 0, 0, 0.1], 1.4),
        t = 0.2, x = [-0.5, 0.5], n = 1001, csv = "csv/sod.csv"))
    println("  wrote csv/{paper,briowu,sod}.csv")

    # 50 random problems: sample(sol, ξ) on 401 points of ξ ∈ [-3, 3]
    P = readdlm(joinpath(DATA, "random_problems.csv"), ',')
    ξ = collect(range(-3, 3; length = 401))
    random = map(1:50) do i
        L, R = P[i, 2:9], P[i, 10:17]
        sol = solve(RiemannProblem(L, R; γ = 5 / 3))
        (seed = Int(P[i, 1]), problem = probrec(L, R, 5 / 3), retcode = sol.retcode,
            speeds = [[r.s_left, r.s_right] for r in wavetable(sol)],
            xi = ξ, W = [collect(sample(sol, z)) for z in ξ])
    end
    save("output", (examples = samples, random = random))
end
