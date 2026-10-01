import { describe, expect, it } from "vitest";
import { juliaRange, juliaRat, juliaRepr } from "../../src/solver/format.ts";
import { ratToDouble } from "../../src/solver/exact.ts";
import { loadGolden } from "../golden.ts";

interface FormatGolden {
  reprs: { x: number; repr: string }[];
  ranges: { a: number; b: number; n: number; x: number[]; repr: string[] }[];
}
const G = loadGolden<FormatGolden>("format");

describe("juliaRepr", () => {
  it("agrees with Julia's repr on all golden values", () => {
    const bad = G.reprs.filter((r) => juliaRepr(r.x) !== r.repr).map((r) => `${r.repr} → ${juliaRepr(r.x)}`);
    expect(bad).toEqual([]);
    expect(G.reprs.length).toBeGreaterThan(500);
  });
  it("handles the non-finite and signed-zero cases", () => {
    expect([NaN, Infinity, -Infinity, 0, -0].map(juliaRepr)).toEqual(["NaN", "Inf", "-Inf", "0.0", "-0.0"]);
  });
});

describe("juliaRange", () => {
  it.each(G.ranges.map((r) => [r.a, r.b, r.n, r] as const))("range(%d, %d; length = %d)", (a, b, n, r) => {
    const x = juliaRange(a, b, n);
    expect(x.length).toBe(n);
    const bad = x.flatMap((v, i) => (Object.is(v, r.x[i]) ? [] : [`${i}: ${v} vs ${r.x[i]}`]));
    expect(bad).toEqual([]);
    expect(x.map(juliaRepr)).toEqual(r.repr);
  });
  it("rejects length 1 with different endpoints, as Julia", () => {
    expect(() => juliaRange(-1, 1, 1)).toThrow(/endpoints differ/);
  });
  it("finds Julia's rational approximations", () => {
    expect(juliaRat(0.1)).toEqual([1, 10]);
    expect(juliaRat(-0.5)).toEqual([1, -2]);
    expect(juliaRat(7.3)).toEqual([73, 10]);
  });
});

describe("ratToDouble", () => {
  it("rounds correctly beyond 2^53", () => {
    expect(ratToDouble(2n ** 60n + 1n, 2n ** 7n)).toBe(2 ** 53);
    expect(ratToDouble(10n ** 30n, 3n * 10n ** 29n)).toBe(10 / 3);
    expect(ratToDouble(-(10n ** 40n), 7n * 10n ** 39n)).toBe(-10 / 7);
    expect(ratToDouble(1n, 10n ** 320n)).toBe(1e-320);
    expect(ratToDouble(3n * 2n ** 100n + 1n, 2n ** 101n)).toBe(1.5);
  });
});
