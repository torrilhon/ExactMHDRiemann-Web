# Shared helpers of the golden-data generator: a minimal JSON writer with exact number
# formatting and serializers for the solver's types.
#
# Numbers are written with `repr` (shortest round-trip, so JSON.parse in JS recovers the
# identical Float64). JSON has no Inf/NaN: they are written as the strings "Inf", "-Inf",
# "NaN" (plan §4.7).

using ExactMHDRiemannSolver, StaticArrays, Pkg, Dates
const E = ExactMHDRiemannSolver

const OUT = get(ENV, "GOLDEN_OUT", normpath(joinpath(@__DIR__, "..", "..", "test", "golden")))
const DATA = joinpath(pkgdir(E), "test", "data")

# ---------------------------------------------------------------- JSON writer

jnum(x::Integer) = string(x)
function jnum(x::Real)
    x = Float64(x)
    isnan(x) && return "\"NaN\""
    isinf(x) && return x > 0 ? "\"Inf\"" : "\"-Inf\""
    return repr(x)
end

function jstr(s::AbstractString)
    io = IOBuffer()
    print(io, '"')
    for c in s
        if c == '"'
            print(io, "\\\"")
        elseif c == '\\'
            print(io, "\\\\")
        elseif c == '\n'
            print(io, "\\n")
        elseif c < ' '
            print(io, "\\u", string(UInt16(c), base = 16, pad = 4))
        else
            print(io, c)
        end
    end
    print(io, '"')
    return String(take!(io))
end

jwrite(io, ::Nothing) = print(io, "null")
jwrite(io, x::Bool) = print(io, x ? "true" : "false")
jwrite(io, x::Real) = print(io, jnum(x))
jwrite(io, x::AbstractString) = print(io, jstr(x))
jwrite(io, x::Symbol) = print(io, jstr(String(x)))
jwrite(io, x::Enum) = print(io, jstr(string(x)))
function jwrite(io, x::Union{AbstractVector,Tuple})
    print(io, '[')
    for (i, v) in enumerate(x)
        i > 1 && print(io, ',')
        jwrite(io, v)
    end
    print(io, ']')
end
function jwrite(io, x::AbstractMatrix)          # row-major list of rows
    jwrite(io, [x[i, :] for i in 1:size(x, 1)])
end
function jwrite_pairs(io, ps)
    print(io, '{')
    for (i, (k, v)) in enumerate(ps)
        i > 1 && print(io, ',')
        print(io, jstr(string(k)), ':')
        jwrite(io, v)
    end
    print(io, '}')
end
jwrite(io, x::NamedTuple) = jwrite_pairs(io, pairs(x))
jwrite(io, x::AbstractDict) = jwrite_pairs(io, sort!(collect(x), by = p -> string(p[1])))

# ---------------------------------------------------------------- provenance

function meta()
    deps = Pkg.dependencies()
    me = deps[Base.UUID("0d315d94-8d33-4e3b-9f90-006d7b2d27bf")]
    names = ("NonlinearSolve", "NonlinearSolveFirstOrder", "NonlinearSolveBase", "OrdinaryDiffEqVerner",
        "OrdinaryDiffEqCore", "Roots", "ForwardDiff", "SciMLBase", "StaticArrays")
    vers = Dict(d.name => string(d.version) for d in values(deps) if d.name in names)
    return (commit = me.git_revision, version = string(me.version), julia = string(VERSION),
        packages = vers, generated = string(Dates.today()))
end

"""Write `data` (a NamedTuple) as test/golden/<name>.json.gz with a `meta` entry first."""
function save(name, data::NamedTuple)
    mkpath(OUT)
    path = joinpath(OUT, name * ".json")
    open(path, "w") do io
        jwrite(io, merge((meta = meta(),), data))
        println(io)
    end
    run(`gzip -9 -n -f $path`)             # deterministic: no name or time stamp in the header
    path *= ".gz"
    println("  wrote ", relpath(path, OUT), " (", round(filesize(path) / 1e6, digits = 2), " MB)")
    return path
end

# ---------------------------------------------------------------- serializers

ser(h::E.HState) = (rho = h.ρ, u = h.u, p = h.p, bt = h.bt, phi = h.φ, vt = [h.vt[1], h.vt[2]])
ser(::Nothing) = nothing
ser(F::E.Frame) = (mirror = F.mirror, sB = F.sB, theta = F.θ, vtR = [F.vtR[1], F.vtR[2]], rho0 = F.ρ0, p0 = F.p0)
ser(c::E.Ctx) = (gamma = c.γ, kappa = c.κ, Bn = c.Bn, scale = collect(c.scale))
ser(f::E.FanData) = (family = f.family, sigma = f.σ, up = ser(f.up), par = f.par)
ser(w::E.Wave) = (kind = w.kind, side = w.side, s_left = w.s_left, s_right = w.s_right,
    left = ser(w.left), right = ser(w.right), fan = ser(w.fan))
sercheck(::Nothing) = nothing
sercheck(c) = (ok = c.ok, maxerr = c.maxerr, messages = c.messages)

"""Wave table row as a compact array: [kind, s_left, s_right, ρ, vx, vy, vz, Bx, By, Bz, p]."""
row(r) = Any[r.kind, r.s_left, r.s_right, r.ρ, r.vx, r.vy, r.vz, r.Bx, r.By, r.Bz, r.p]

"""Full record of a solve. `full = false` drops the canonical waves (only the user-frame table)."""
function solrec(sol::RiemannSolution; full::Bool = true, time_ms = nothing)
    has = !isempty(sol.waves) && sol.frame !== nothing
    base = (retcode = sol.retcode, reason = sol.reason, method = sol.method, psi = collect(sol.Ψ),
        residual = sol.residual, check = sercheck(sol.check),
        table = has ? [row(r) for r in wavetable(sol)] : Any[], time_ms = time_ms)
    full || return base
    return merge(base, (frame = ser(sol.frame), ctx = sol.ctx === nothing ? nothing : ser(sol.ctx),
        waves = [ser(w) for w in sol.waves]))
end

optsrec(o::SolverOptions) = Dict(string(k) => getfield(o, k) for k in fieldnames(SolverOptions))

"""Solve with `time_limit = Inf` (plan §4.2) and return (record, sol)."""
function golden_solve(L, R, γ; opts = SolverOptions(), full = true)
    o = SolverOptions(; (k => getfield(opts, k) for k in fieldnames(SolverOptions))..., time_limit = Inf)
    t = @elapsed sol = solve(RiemannProblem(L, R; γ), opts = o)
    return solrec(sol; full, time_ms = 1e3t), sol
end

probrec(L, R, γ) = (L = collect(Float64, L), R = collect(Float64, R), gamma = Float64(γ))
