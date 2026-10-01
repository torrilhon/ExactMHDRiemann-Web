// L2 golden tests (plan §5.1, §5.3): the 5×5 residual at 1400 Ψ (including domain
// edges and branch switches) and the finite-difference Jacobian against ForwardDiff.

import { describe, expect, it } from "vitest";
import { isBig, residualJacobian, safeResidual } from "../src/solver/residual.ts";
import { defaultOptions, type Ctx } from "../src/solver/types.ts";
import { loadGolden } from "./golden.ts";
import { fromJson, relerr, Worst } from "./helpers.ts";

const G = loadGolden<any>("residual");

function ctxOf(c: any): Ctx {
  return { gamma: c.gamma, kappa: c.kappa, Bn: c.Bn, opts: defaultOptions(), scale: c.scale };
}

describe("L2 residual", () => {
  it("agrees with Julia at all points; BIG exactly where Julia returns it", () => {
    const w = new Worst(1e-10);
    const bigMismatch: string[] = [];
    let nbig = 0;
    G.problems.forEach((p: any, i: number) => {
      const UL = fromJson(p.UL), UR = fromJson(p.UR), ctx = ctxOf(p.ctx);
      p.evals.forEach((e: any, j: number) => {
        const r = safeResidual(e.psi, UL, UR, ctx);
        if (isBig(r) !== e.big) bigMismatch.push(`${i}.${j} Ψ=${JSON.stringify(e.psi)}: ts ${isBig(r)} julia ${e.big}`);
        if (e.big) nbig++;
        else w.add(`${i}.${j}`, relerr(r, e.r));
      });
    });
    expect(bigMismatch).toEqual([]);
    expect(nbig).toBe(130);
    expect(w.bad).toEqual([]);
    console.log(`residual: worst relerr ${w.max.toExponential(2)} over ${G.problems.length * 56 - nbig} points`);
  });

  it("finite-difference Jacobian agrees with ForwardDiff to 1e-6 (runtests.jl criterion)", () => {
    const w = new Worst(1e-6);
    G.problems.forEach((p: any, i: number) => {
      const UL = fromJson(p.UL), UR = fromJson(p.UR), ctx = ctxOf(p.ctx);
      p.jacobians.forEach((jc: any, j: number) => {
        const J = residualJacobian(jc.psi, UL, UR, ctx);
        const scale = Math.max(1, ...jc.J.flat().map(Math.abs));
        let e = 0;
        J.forEach((row, a) => row.forEach((v, b) => { e = Math.max(e, Math.abs(v - jc.J[a][b])); }));
        w.add(`${i}.${j} Ψ=${JSON.stringify(jc.psi.map((x: number) => +x.toPrecision(3)))}`, e / scale);
      });
    });
    console.log(`Jacobian: worst ${w.max.toExponential(2)}; failing: ${w.bad.length}`);
    expect(w.bad).toEqual([]);
  });
});
