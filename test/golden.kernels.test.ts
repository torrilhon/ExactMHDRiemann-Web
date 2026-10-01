// L0 and L1 golden tests (plan §5.1, §5.3): basics, canonical frame, shock and fan
// kernels, dense fan output, against the Julia reference.

import { describe, expect, it } from "vitest";
import { conserved, flux, primitiveJacobian, speedsW } from "../src/solver/eos.ts";
import { canonicalStates, fromCanonical, makeFrame, speedToUser, toHState } from "../src/solver/canonical.ts";
import { fastCubicBt, fastShock, largestRootIn, slowDefect, slowLimit, slowShockState, slowVolume } from "../src/solver/shocks.ts";
import { fanSpeed, fanState, fastFan, slowFan } from "../src/solver/fans.ts";
import { rotation } from "../src/solver/side.ts";
import { defaultOptions, makeCtx } from "../src/solver/types.ts";
import { loadGolden } from "./golden.ts";
import { fromJson, hrelerr, relerr, Worst, type HStateJson } from "./helpers.ts";

const flat = (m: number[][]) => m.flat();

describe("L0 basics", () => {
  const G = loadGolden<{ states: any[]; frames: any[] }>("basics");
  it("speeds, conserved, flux, primitive Jacobian (≤ 1e-14)", () => {
    const w = new Worst(1e-14);
    G.states.forEach((s, i) => {
      w.add(`speeds ${i}`, relerr(speedsW(s.W, s.gamma), s.speeds));
      w.add(`conserved ${i}`, relerr(conserved(s.W, s.gamma), s.conserved));
      w.add(`flux ${i}`, relerr(flux(s.W, s.gamma), s.flux));
      w.add(`jacobian ${i}`, relerr(flat(primitiveJacobian(s.W, s.gamma)), flat(s.jacobian)));
    });
    expect(w.bad).toEqual([]);
  });
  it("canonical frame and its inverse (≤ 1e-14)", () => {
    const w = new Worst(1e-14);
    G.frames.forEach((f, i) => {
      const F = makeFrame(f.L, f.R);
      expect(F.mirror).toBe(f.frame.mirror);
      expect(F.sB).toBe(f.frame.sB);
      w.add(`frame ${i}`, relerr([F.theta, ...F.vtR, F.rho0, F.p0], [f.frame.theta, ...f.frame.vtR, f.frame.rho0, f.frame.p0]));
      const [Lc, Rc] = canonicalStates(f.L, f.R, F);
      w.add(`Lc ${i}`, relerr(Lc, f.Lc));
      w.add(`Rc ${i}`, relerr(Rc, f.Rc));
      w.add(`UL ${i}`, hrelerr(toHState(Lc), f.UL));
      w.add(`UR ${i}`, hrelerr(toHState(Rc), f.UR));
      const Lb = F.mirror ? fromCanonical(Rc, F) : fromCanonical(Lc, F);
      const Rb = F.mirror ? fromCanonical(Lc, F) : fromCanonical(Rc, F);
      w.add(`back ${i}`, relerr([...Lb, ...Rb], [...f.Lback, ...f.Rback]));
      w.add(`speed ${i}`, relerr([-1.3, 0, 0.7].map((s) => speedToUser(s, F)), f.speed_to_user));
    });
    expect(w.bad).toEqual([]);
  });
});

describe("L1 kernels", () => {
  const G = loadGolden<any>("kernels");
  const ctxOf = (c: { gamma: number; Bn: number }, o = {}) => makeCtx(c.gamma, c.Bn, defaultOptions(o));

  it("fast shocks: cubic roots (≥ 95 % bit-identical, else ≤ 8 ULP), states ≤ 1e-13", () => {
    const w = new Worst(1e-13);
    const roots: string[] = [];
    const all = [...G.shocks.flatMap((s: any) => s.fast.map((f: any) => ({ ctx: s.ctx, U: s.U, f, o: {} }))),
      ...G.switch_on.map((s: any) => ({ ctx: s.ctx, U: s.U, f: s.fast, o: {} })),
      ...G.small_bn.map((s: any) => ({ ctx: s.ctx, U: s.U, f: s.fast, o: { bn_euler: 0 } }))];
    all.forEach(({ ctx, U, f, o }: any, i: number) => {
      const c = ctxOf(ctx, o);
      const cub = f.cubic;
      // the coefficients as the solver computes them, and the root of Brent
      const r = largestRootIn(cub.c, cub.lo, cub.hi);
      if (!Object.is(r, cub.root)) {
        roots.push(`${i}: ${r} vs ${cub.root}`);
        // Julia's evalpoly is partly fused by LLVM: the root may differ in the last bits
        // (Brent stops at 4 ULP)
        w.add(`root ${i}`, Math.abs(r - cub.root) / Math.abs(cub.root) <= 8 * Number.EPSILON ? 0 : 1);
      }
      const [D, s] = fastShock(fromJson(U), f.psi, f.sigma, c);
      w.add(`fast ${i} ψ=${f.psi}`, hrelerr(D, f.down, [s], [f.speed]));
    });
    expect(roots.length).toBeLessThanOrEqual(0.05 * all.length);
    expect(w.bad).toEqual([]);
    expect(all.length).toBe(120 * 5 + 12 + 9);
    console.log(`cubic roots: ${all.length - roots.length}/${all.length} bit-identical; states worst ${w.max.toExponential(2)}`);
  });

  it("slow limit, slow shocks, rotations (≤ 1e-13)", () => {
    const w = new Worst(1e-13);
    G.shocks.forEach((s: any, i: number) => {
      const c = ctxOf(s.ctx);
      const U = fromJson(s.U);
      w.add(`slow_limit ${i}`, relerr([slowLimit(U, c)], [s.slow_limit]));
      s.slow.forEach((sl: any, j: number) => {
        const [D, sp, M, v] = slowShockState(U, sl.delta, s.sigma, c);
        w.add(`slow ${i}.${j}`, hrelerr(D, sl.down, [sp, M, v], [sl.speed, sl.M, sl.v]));
        const A = U.bt / Math.sqrt(U.p), B2 = (c.Bn / Math.sqrt(U.p)) ** 2;
        w.add(`volume ${i}.${j}`, relerr(slowVolume(sl.delta, A, B2, c.kappa), sl.volume));
        w.add(`defect ${i}.${j}`, relerr([slowDefect(U, sl.delta, c)], [sl.defect]));
      });
      const [Dr, sr] = rotation(U, s.rotation.alpha, s.sigma, c);
      w.add(`rotation ${i}`, hrelerr(Dr, s.rotation.down, [sr], [s.rotation.speed]));
    });
    G.slow_small.forEach((s: any, i: number) => {
      w.add(`slow_limit small ${i}`, relerr([slowLimit(fromJson(s.U), ctxOf(s.ctx))], [s.slow_limit]));
    });
    expect(w.bad).toEqual([]);
  });

  it("fan end states (≤ 1e-10)", () => {
    const w = new Worst(1e-10);
    G.fans.forEach((f: any, i: number) => {
      const c = ctxOf(f.ctx);
      const U = fromJson(f.U);
      w.add(`fast fan ${i}`, hrelerr(fastFan(U, f.fast.psi, f.sigma, c)[0], f.fast.down));
      w.add(`slow fan ${i}`, hrelerr(slowFan(U, f.slow.psi, f.sigma, c)[0], f.slow.down));
    });
    expect(w.bad).toEqual([]);
    console.log(`fan end states: worst relerr ${w.max.toExponential(2)}`);
  });
});

describe("L1 dense fans", () => {
  const G = loadGolden<{ fans: any[] }>("fans_dense");
  it("fan_state and fan_speed at 50 τ per fan (≤ 1e-10)", () => {
    const w = new Worst(1e-10);
    let n = 0;
    for (const [i, f] of G.fans.entries()) {
      if (f.family === "fast0") continue;                       // quasi-Euler: Phase 7
      const c = makeCtx(f.ctx.gamma, f.ctx.Bn, defaultOptions());
      const U = fromJson(f.up);
      const [, fan] = f.family === "fast" ? fastFan(U, f.psi, f.sigma, c, true) : slowFan(U, f.psi, f.sigma, c, true);
      f.tau.forEach((tau: number, k: number) => {
        w.add(`${f.family} ${i} τ=${tau}`, hrelerr(fanState(fan!, tau, c), f.states[k] as HStateJson,
          [fanSpeed(fan!, tau, c)], [f.speeds[k]]));
        n++;
      });
    }
    expect(w.bad).toEqual([]);
    expect(n).toBe(100 * 50);
    console.log(`dense fans: worst relerr ${w.max.toExponential(2)}`);
  });
});
