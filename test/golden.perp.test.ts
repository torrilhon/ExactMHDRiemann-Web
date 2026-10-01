// Quasi-Euler kernels (perpendicular.jl) against the Julia reference.

import { describe, expect, it } from "vitest";
import { fanSpeed, fanState } from "../src/solver/fans.ts";
import { perpFan, perpShock, perpSide } from "../src/solver/perpendicular.ts";
import { DomainError } from "../src/solver/math.ts";
import { defaultOptions, makeCtx } from "../src/solver/types.ts";
import { loadGolden } from "./golden.ts";
import { fromJson, hrelerr, relerr, Worst, type HStateJson } from "./helpers.ts";

describe("quasi-Euler kernels", () => {
  const K = loadGolden<any>("kernels");
  it("perp_side for P above and below the total pressure (≤ 1e-12)", () => {
    const w = new Worst(1e-12);
    K.perp.forEach((c: any, i: number) => {
      const ctx = makeCtx(c.ctx.gamma, c.ctx.Bn, defaultOptions());
      const [D, wave] = perpSide(fromJson(c.U), c.P, c.sigma, ctx, true);
      w.add(`perp ${i}`, hrelerr(D, c.down));
      expect(wave?.kind ?? null).toBe(c.wave?.kind ?? null);
      if (wave && c.wave) {
        w.add(`perp speeds ${i}`, relerr([wave.s_left, wave.s_right], [c.wave.s_left, c.wave.s_right]));
        if (c.wave.send !== null) w.add(`perp send ${i}`, relerr([wave.fan!.par], [c.wave.send]));
      }
    });
    expect(w.bad).toEqual([]);
    expect(K.perp.length).toBe(660);
  });
  it("refuses shocks beyond P/P0 ~ 1e14 as Julia does", () => {
    for (const c of K.perp_refusal) {
      const ctx = makeCtx(c.ctx.gamma, c.ctx.Bn, defaultOptions());
      let got: string;
      try { perpShock(fromJson(c.U), c.P, 1, ctx); got = "ok"; } catch (e) {
        if (!(e instanceof DomainError)) throw e;
        got = "DomainError";
      }
      expect(got).toBe(typeof c.result === "string" ? c.result : "ok");
    }
  });
  it("dense quasi-Euler fans (:fast0) at 50 τ (≤ 1e-10)", () => {
    const D = loadGolden<{ fans: any[] }>("fans_dense");
    const w = new Worst(1e-10);
    let n = 0;
    for (const [i, f] of D.fans.entries()) {
      if (f.family !== "fast0") continue;
      const ctx = makeCtx(f.ctx.gamma, f.ctx.Bn, defaultOptions());
      const [, fan] = perpFan(fromJson(f.up), f.P, f.sigma, ctx, true);
      f.tau.forEach((tau: number, k: number) => {
        w.add(`fast0 ${i} τ=${tau}`, hrelerr(fanState(fan!, tau, ctx), f.states[k] as HStateJson, [fanSpeed(fan!, tau, ctx)], [f.speeds[k]]));
        n++;
      });
    }
    expect(n).toBe(20 * 50);
    expect(w.bad).toEqual([]);
  });
});

import { wavetable } from "../src/solver/output.ts";
import { solve } from "../src/solver/solve.ts";
import { riemannProblem } from "../src/solver/types.ts";

describe("Sod (Toro, Riemann Solvers, Table 4.2 test 1)", () => {
  it("matches Toro's values to 1e-5", () => {
    const sol = solve(riemannProblem([1, 0, 0, 0, 0, 0, 0, 1], [0.125, 0, 0, 0, 0, 0, 0, 0.1], 1.4));
    expect(sol.retcode).toBe("Success");
    expect(sol.method).toBe("quasi_euler");
    const T = wavetable(sol);
    expect(T.map((r) => r.kind)).toEqual(["fast_fan", "tangential", "fast_shock"]);
    expect(Math.abs(T[1]!.p - 0.30313)).toBeLessThan(1e-5);
    expect(Math.abs(T[1]!.vx - 0.92745)).toBeLessThan(1e-5);
    expect(Math.abs(T[0]!.rho - 0.42632)).toBeLessThan(1e-5);
    expect(Math.abs(T[1]!.rho - 0.26557)).toBeLessThan(1e-5);
    expect(Math.abs(T[2]!.s_left - 1.75216)).toBeLessThan(1e-5);
  });
});
