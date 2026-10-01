# L0: characteristic speeds, conserved variables, fluxes, primitive Jacobian, canonical
# frame; L0 format: Julia's number printing and ranges (plan §4.7).

function randW(rng; Bx = 4rand(rng) - 2)
    [0.1 + 5rand(rng), 4rand(rng) - 2, 4rand(rng) - 2, 4rand(rng) - 2, Bx, 4rand(rng) - 2, 4rand(rng) - 2, 0.1 + 5rand(rng)]
end

function gen_basics()
    rng = MersenneTwister(11)
    states = map(1:500) do _
        W = randW(rng)
        rand(rng) < 0.1 && (W[6] = W[7] = 0.0)          # Bt = 0
        rand(rng) < 0.1 && (W[5] = 0.0)                 # Bn = 0
        γ = rand(rng, (1.4, 5 / 3, 2.0))
        (W = W, gamma = γ, speeds = collect(E.speeds(SVector{8}(W), γ)),
            conserved = collect(E.conserved(W, γ)), flux = collect(E.flux(W, γ)),
            jacobian = E.primitive_jacobian(W, γ))
    end
    # canonical frame: random pairs with equal Bx, incl. mirrored and Bx < 0 cases
    frames = map(1:300) do i
        Bx = (2rand(rng) - 1) * 2
        L, R = randW(rng; Bx), randW(rng; Bx)
        i % 10 == 0 && (R[6:7] = -0.7 .* L[6:7])                 # coplanar, antiparallel
        F = E.make_frame(SVector{8}(L), SVector{8}(R))
        Lc, Rc = E.canonical_states(SVector{8}(L), SVector{8}(R), F)
        Lb = F.mirror ? E.from_canonical(Rc, F) : E.from_canonical(Lc, F)
        Rb = F.mirror ? E.from_canonical(Lc, F) : E.from_canonical(Rc, F)
        (L = L, R = R, frame = ser(F), Lc = collect(Lc), Rc = collect(Rc),
            UL = ser(E.to_hstate(Lc)), UR = ser(E.to_hstate(Rc)), Lback = collect(Lb), Rback = collect(Rb),
            speed_to_user = [E.speed_to_user(s, F) for s in (-1.3, 0.0, 0.7)])
    end
    save("basics", (states = states, frames = frames))
end

function gen_format()
    rng = MersenneTwister(12)
    special = [0.0, -0.0, 1.0, -1.0, 0.1, 0.5, 1 / 3, 2 / 3, 1.5, 100.0, 1e5, 1e6, 123456.789,
        1e15, 1e16, 1.0e17, 1.0e21, 1.0e22, 9007199254740993.0, 12345678901234567.0,
        1e-3, 1e-4, 1e-5, 0.000123, 1.5e-5, 1e-7, 1e-300, 5e-324, 2.2250738585072014e-308,
        floatmax(), prevfloat(1.0), nextfloat(1.0), eps(), π, -π, 2e-10, 0.9e-10, 1 / 7, 1e100, -2.5e-8]
    for _ in 1:400
        push!(special, (rand(rng) < 0.5 ? -1 : 1) * rand(rng) * 10.0^rand(rng, -20:20))
    end
    for e in -8:24                       # around the switch between fixed and scientific notation
        push!(special, 10.0^e, 1.234 * 10.0^e, prevfloat(10.0^e), nextfloat(10.0^e))
    end
    reprs = [(x = x, repr = repr(x)) for x in special]
    grids = [(-1.0, 1.0, 2001), (-0.5, 0.5, 2001), (-1.0, 1.0, 101), (-1.0, 3.0, 7), (-1.0, 3.0, 1000),
        (0.0, 1.0, 3), (-2.5, 7.3, 333), (0.1, 0.3, 11), (-3.0, 3.0, 61), (1e-3, 2e-3, 17), (5.0, -5.0, 9),
        (2.0, 2.0, 4), (-1.0, -1.0, 1)]
    ranges = [(a = a, b = b, n = n, x = collect(range(a, b; length = n)), repr = [repr(v) for v in range(a, b; length = n)])
              for (a, b, n) in grids]
    save("format", (reprs = reprs, ranges = ranges))
end
