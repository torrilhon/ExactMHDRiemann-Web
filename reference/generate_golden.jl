# Generates the golden data in test/golden/ (plan §5.1). Run from the repository root:
#
#   julia --project=reference reference/generate_golden.jl            # everything
#   julia --project=reference reference/generate_golden.jl kernels    # selected parts
#
# Parts: basics format kernels fans_dense residual solve_examples solve_random
#        solve_stress solve_small_bt solve_perp output data

using Random
include("golden/common.jl")
include("golden/basics.jl")
include("golden/kernels.jl")
include("golden/residual.jl")
include("golden/solves.jl")
include("golden/scans.jl")
include("golden/output.jl")

function copy_data()
    mkpath(joinpath(OUT, "data"))
    for f in ("random_problems.csv", "c_reference.csv", "README.md")
        cp(joinpath(DATA, f), joinpath(OUT, "data", f); force = true)
    end
    println("  copied test/data of the Julia package to data/")
end

const PARTS = [
    "basics" => gen_basics, "format" => gen_format, "kernels" => gen_kernels, "fans_dense" => gen_fans_dense,
    "residual" => gen_residual, "solve_examples" => gen_solve_examples, "solve_random" => gen_solve_random,
    "solve_stress" => gen_solve_stress, "solve_small_bt" => gen_small_bt_scan, "solve_perp" => gen_solve_perp,
    "output" => gen_output, "data" => copy_data]

selected = isempty(ARGS) ? first.(PARTS) : ARGS
for a in selected
    any(p -> p[1] == a, PARTS) || error("unknown part $a; parts: $(join(first.(PARTS), ' '))")
end
println("ExactMHDRiemannSolver ", meta().commit, ", Julia ", VERSION)
for (name, f) in PARTS
    name in selected || continue
    println(name, ":")
    t = @elapsed f()
    println("  ", round(t, digits = 1), " s")
end
