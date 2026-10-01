import { describe, expect, it } from "vitest";
import { fdJacobian } from "../../src/solver/jacobian.ts";

describe("finite-difference Jacobian (plan §4.3)", () => {
  const f = (x: readonly number[]) => [Math.sin(x[0]!) * x[1]!, x[0]! ** 3 + Math.exp(x[1]!)];
  it("is accurate to ~1e-10 with central differences", () => {
    const x = [0.7, -0.4];
    const { J, ok } = fdJacobian(f, x, f(x));
    const ex = [[Math.cos(0.7) * -0.4, Math.sin(0.7)], [3 * 0.49, Math.exp(-0.4)]];
    expect(ok).toEqual([true, true]);
    J.forEach((r, i) => r.forEach((v, j) => expect(Math.abs(v - ex[i]![j]!)).toBeLessThan(1e-9)));
  });
  it("stays on the base point's branch at a kink", () => {
    // g has slope 1 for x ≤ 0 and slope 3 for x > 0
    const g = (x: readonly number[]) => [x[0]! <= 0 ? x[0]! : 3 * x[0]!];
    const at = (x0: number) => fdJacobian(g, [x0], g([x0]), { switches: [[0]] }).J[0]![0]!;
    expect(at(0)).toBeCloseTo(1, 7);
    expect(at(-1e-8)).toBeCloseTo(1, 7);
    expect(at(1e-8)).toBeCloseTo(3, 7);
    expect(at(0.5)).toBeCloseTo(3, 9);
    // without the switch information the stencil averages the slopes
    expect(fdJacobian(g, [0], [0]).J[0]![0]).toBeCloseTo(2, 7);
  });
  it("uses the other side at a domain edge and marks dead columns", () => {
    const BIG = 1e6;
    // valid for x0 ≤ 1; in x1 valid only exactly at 6 (both perturbations invalid)
    const h = (x: readonly number[]) => (x[0]! > 1 || x[1]! !== 6 ? [BIG, BIG] : [x[0]! * 2, x[1]! * x[1]!]);
    const inv = (r: readonly number[]) => r.every((v) => v === BIG);
    const x = [1, 6];
    const { J, ok } = fdJacobian(h, x, h(x), { invalid: inv });
    expect(ok).toEqual([true, false]);
    expect(J[0]![0]).toBeCloseTo(2, 9);              // backward difference
    expect(J[1]![1]).toBe(0);
  });
  it("calls beforePerturbed once per perturbed evaluation", () => {
    let calls = 0, evals = 0;
    const k = (x: readonly number[]) => { evals++; return [x[0]!, x[1]!]; };
    fdJacobian(k, [1, 2], [1, 2], { beforePerturbed: () => calls++ });
    expect(calls).toBe(evals);
    expect(evals).toBe(4);
  });
});
