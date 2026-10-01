// L4 golden tests (plan §5.3): sampled profiles and the CSV export.

import { describe, expect, it } from "vitest";
import { sampleXi, toCSV } from "../src/solver/output.ts";
import { solve } from "../src/solver/solve.ts";
import { riemannProblem } from "../src/solver/types.ts";
import { goldenText, loadGolden } from "./golden.ts";
import { relerr, Worst } from "./helpers.ts";

const G = loadGolden<any>("output");

describe("L4 sampled profiles", () => {
  it("agree with Julia to 1e-10 away from the wave speeds (50 random problems × 401 ξ)", () => {
    const w = new Worst(1e-10);
    let skipped = 0;
    for (const c of G.random) {
      const sol = solve(riemannProblem(c.problem.L, c.problem.R, c.problem.gamma));
      expect(sol.retcode).toBe(c.retcode);
      const edges: number[] = c.speeds.flat();
      c.xi.forEach((xi: number, k: number) => {
        // points within 1e-12 of a wave speed may fall on different sides
        if (edges.some((s) => Math.abs(xi - s) <= 1e-12 * Math.max(1, Math.abs(s)))) { skipped++; return; }
        w.add(`seed ${c.seed} ξ=${xi}`, relerr(sampleXi(sol, xi), c.W[k]));
      });
    }
    expect(w.bad).toEqual([]);
    console.log(`profiles: worst relerr ${w.max.toExponential(2)}, skipped ${skipped} points at wave speeds`);
  });
});

describe("CSV export (write_csv)", () => {
  it.each(G.examples as any[])("$name: same layout, x column identical, values ≤ 1e-10", (e: any) => {
    const sol = solve(riemannProblem(e.problem.L, e.problem.R, e.problem.gamma));
    const ts = toCSV(sol, e.t, e.x[0], e.x[1], e.n).split("\n");
    const jl = goldenText(e.csv).split("\n");
    expect(ts.length).toBe(jl.length);
    expect(ts[0]).toBe(jl[0]);
    expect(ts[ts.length - 1]).toBe("");                          // trailing newline as println
    const w = new Worst(1e-10);
    let identical = 0;
    for (let i = 1; i < jl.length - 1; i++) {
      const a = ts[i]!.split(","), b = jl[i]!.split(",");
      expect(a[0]).toBe(b[0]);                                   // juliaRange + juliaRepr
      if (ts[i] === jl[i]) identical++;
      w.add(`line ${i}`, relerr(a.slice(1).map(Number), b.slice(1).map(Number)));
    }
    expect(w.bad).toEqual([]);
    console.log(`${e.name}.csv: ${identical}/${jl.length - 2} lines character-identical, worst relerr ${w.max.toExponential(2)}`);
  });
});
