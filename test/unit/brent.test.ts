import { describe, expect, it } from "vitest";
import { ConvergenceFailed, findZeroBrent, nextfloat, prevfloat } from "../../src/solver/brent.ts";
import { ArgumentError, EPS } from "../../src/solver/math.ts";
import { evalpolyFma } from "../fma.ts";
import { loadGolden } from "../golden.ts";

interface BrentCase {
  name: string; c: number[]; a: number; b: number;
  opts: { xatol?: number; xrtol?: number }; root: number | string; evals: number;
}
const G = loadGolden<{ cases: BrentCase[] }>("brent");

function run(c: BrentCase): { root: number | string; evals: number } {
  let evals = 0;
  const f = (x: number) => { evals++; return evalpolyFma(x, c.c); };       // golden made with fused evalpoly
  try {
    return { root: findZeroBrent(f, c.a, c.b, c.opts), evals };
  } catch (e) {
    if (e instanceof ConvergenceFailed) return { root: "ConvergenceFailed", evals };
    if (e instanceof ArgumentError) return { root: "ArgumentError", evals };
    throw e;
  }
}

describe("Brent vs Roots.jl (bit level, plan §4.4)", () => {
  it("has golden cases", () => expect(G.cases.length).toBeGreaterThan(300));
  it("reproduces every root and evaluation count exactly", () => {
    const bad = G.cases.flatMap((c) => {
      const r = run(c);
      return Object.is(r.root, c.root) && r.evals === c.evals
        ? [] : [`${c.name}: ${String(r.root)} (${r.evals}) vs ${String(c.root)} (${c.evals})`];
    });
    expect(bad).toEqual([]);
  });
});

describe("Brent edge cases", () => {
  it("rejects a bracket without sign change", () => {
    expect(() => findZeroBrent((x) => x * x + 1, -1, 1)).toThrow(ArgumentError);
  });
  it("reaches 4 ULP on a wide bracket", () => {
    const r = findZeroBrent((x) => x - 1e-8, 1e-300, 1, { xatol: 0, xrtol: 4 * EPS });
    expect(Math.abs(r - 1e-8)).toBeLessThanOrEqual(8 * EPS * 1e-8);
  });
  it("nextfloat / prevfloat", () => {
    expect(nextfloat(1)).toBe(1 + EPS);
    expect(prevfloat(1)).toBe(1 - EPS / 2);
    expect(nextfloat(0)).toBe(Number.MIN_VALUE);
    expect(nextfloat(-Number.MIN_VALUE)).toBe(-0);
  });
});
