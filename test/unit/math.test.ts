import { describe, expect, it } from "vitest";
import { DomainError, evalpoly, logD, powD, sqrtD } from "../../src/solver/math.ts";

describe("checked functions (plan §4.6)", () => {
  it("throw DomainError where Julia does", () => {
    expect(() => sqrtD(-1e-300)).toThrow(DomainError);
    expect(() => logD(-1)).toThrow(DomainError);
    expect(() => powD(-2, 0.5)).toThrow(DomainError);
  });
  it("agree with Julia at the domain boundary", () => {
    expect(Object.is(sqrtD(-0), -0)).toBe(true);    // Julia: sqrt(-0.0) == -0.0
    expect(logD(0)).toBe(-Infinity);
    expect(Number.isNaN(sqrtD(NaN))).toBe(true);
    expect(powD(-2, 3)).toBe(-8);
    expect(powD(0, 0.5)).toBe(0);
  });
  it("evalpoly uses ascending coefficients (Horner)", () => {
    expect(evalpoly(2, [1, 2, 3])).toBe(1 + 4 + 12);
    expect(evalpoly(0.5, [4])).toBe(4);
  });
});

import { fma } from "../../src/solver/math.ts";
import { loadGolden } from "../golden.ts";

describe("fma emulation vs Julia's hardware fma", () => {
  const G = loadGolden<{ fma: { a: number; b: number; c: number; fma: number }[];
    evalpoly: { c: number[]; x: number; value: number }[] }>("format");
  it("is bit-identical on all golden triples", () => {
    const bad = G.fma.filter((t) => !Object.is(fma(t.a, t.b, t.c), t.fma))
      .map((t) => `fma(${t.a}, ${t.b}, ${t.c}) = ${fma(t.a, t.b, t.c)}, Julia ${t.fma}`);
    expect(bad).toEqual([]);
    expect(G.fma.length).toBeGreaterThan(3000);
  });
  it("evalpoly is bit-identical to Julia's", () => {
    const bad = G.evalpoly.filter((p) => !Object.is(evalpoly(p.x, p.c), p.value));
    expect(bad.length).toBe(0);
  });
  it("differs from a separately rounded multiply-add where it must", () => {
    // the case found in the Brent comparison: plain Horner gives exactly 0
    const x = 0.7235214403368077, c0 = 0.3074344683174547, c1 = -0.4249141092133225;
    expect(c0 + x * c1).toBe(0);
    expect(evalpoly(x, [c0, c1])).toBe(-3.9650239443466465e-18);
  });
});
