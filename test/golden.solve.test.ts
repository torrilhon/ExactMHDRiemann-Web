// L3 golden tests (plan §5.2, §5.3) for the sets that are fast enough for every CI run:
// examples and edge cases, the 2000 random problems with the C reference, and the
// quasi-Euler set. The stress and small-Bt sets run with scripts/compare-golden.mts.

import { describe, expect, it } from "vitest";
import { sampleXi, wavetable } from "../src/solver/output.ts";
import { caseLabel, runCase, type GoldenCase } from "./compare.ts";
import { goldenText, loadGolden, type GoldenName } from "./golden.ts";

// known, explained differences (test/golden/DIFFERENCES.md): Julia's checker false
// positive on vanishingly weak fans, fixed in the port (Julia commit 0ca567a)
const EXPECTED: Record<string, { retcode: string; reason: string }> = {
  "refusal/bt_larger_at_threshold": { retcode: "Success", reason: "none" },
  "refusal/bt_ok_2": { retcode: "Success", reason: "none" },
};

function runSet(set: GoldenName) {
  const G = loadGolden<{ cases?: GoldenCase[]; scan?: GoldenCase[] }>(set);
  const cases = [...(G.cases ?? []), ...(G.scan ?? [])];
  const codeBad: string[] = [], tableBad: string[] = [], checkBad: string[] = [];
  let worst = 0;
  cases.forEach((c, i) => {
    const r = runCase(c, caseLabel(c, i));
    const exp = EXPECTED[r.label];
    if (exp) {
      if (r.sol.retcode !== exp.retcode || r.sol.reason !== exp.reason) codeBad.push(`${r.label}: ${r.sol.retcode}`);
      return;
    }
    if (!r.codeMatch) codeBad.push(`${r.label}: julia ${c.retcode}/${c.reason}/${c.method} ts ${r.sol.retcode}/${r.sol.reason}/${r.sol.method}`);
    else if (!(r.tableErr <= r.tableTol)) tableBad.push(`${r.label}: ${r.tableErr}`);
    else worst = Math.max(worst, r.tableErr);
    if (r.sol.retcode === "CheckFailed" || (r.sol.retcode === "Success" && !r.checkOk)) checkBad.push(r.label);
  });
  return { n: cases.length, codeBad, tableBad, checkBad, worst };
}

describe.each(["solve_examples", "solve_random", "solve_perp"] as const)("L3 %s", (set) => {
  const res = runSet(set);
  it("identical retcode, reason and method", () => expect(res.codeBad).toEqual([]));
  it("wave tables within max(1e-9, 1e4 × residual)", () => {
    expect(res.tableBad).toEqual([]);
    console.log(`${set}: ${res.n} problems, worst table relerr ${res.worst.toExponential(2)}`);
  });
  it("no CheckFailed, no Success with check messages", () => expect(res.checkBad).toEqual([]));
});

describe("C reference (test/data/c_reference.csv, 8 decimals)", () => {
  // The README's "1.3e-8" is rounded: the deviation is 1.334e-8, the same for Julia and
  // the port (which agree to 8e-15 on these problems); the C data carry 8 decimals.
  it("states left and right of the contact within 1.4e-8 (Julia: 1.3e-8, rounded)", () => {
    const G = loadGolden<{ cases: GoldenCase[] }>("solve_random");
    const bySeed = new Map(G.cases.map((c) => [c.seed!, c]));
    const rows = goldenText("data/c_reference.csv").trim().split("\n").map((l) => l.split(",").map(Number));
    let worst = 0;
    for (const row of rows) {
      const c = bySeed.get(row[0]!)!;
      const r = runCase(c, `seed ${row[0]}`);
      const sc = wavetable(r.sol).find((w) => w.kind === "contact")!.s_left;
      const a = sampleXi(r.sol, sc - 1e-9), b = sampleXi(r.sol, sc + 1e-9);
      a.forEach((v, j) => { worst = Math.max(worst, Math.abs(v - row[17 + j]!)); });
      b.forEach((v, j) => { worst = Math.max(worst, Math.abs(v - row[25 + j]!)); });
    }
    expect(rows.length).toBe(355);
    expect(worst).toBeLessThan(1.4e-8);
    console.log(`C reference: worst ${worst.toExponential(2)}`);
  });
});
