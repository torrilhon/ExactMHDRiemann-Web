import { describe, expect, it } from "vitest";
import { luFactor, luSolve, matTvec, matvec, normInf, opnormInf, solve } from "../../src/solver/linalg.ts";

describe("linalg", () => {
  it("solves a 5×5 system that needs pivoting", () => {
    const A = [
      [0, 2, 1, 0, 3], [1, 0, 4, 1, 0], [2, 1, 0, 5, 1], [0, 3, 1, 1, 0], [1, 1, 1, 1, 1e-3],
    ];
    const x = [1, -2, 0.5, 3, -1];
    const b = matvec(A, x);
    const y = solve(A, b)!;
    y.forEach((v, i) => expect(v).toBeCloseTo(x[i]!, 13));
  });
  it("reports singular matrices", () => {
    expect(solve([[1, 2], [2, 4]], [1, 2])).toBeNull();
    expect(luFactor([[0, 0], [0, 1]]).singular).toBe(true);
  });
  it("is exact on a permutation", () => {
    const f = luFactor([[0, 1], [1, 0]]);
    expect(luSolve(f, [3, 4])).toEqual([4, 3]);
  });
  it("norms", () => {
    expect(normInf([1, -3, 2])).toBe(3);
    expect(Number.isNaN(normInf([1, NaN]))).toBe(true);
    expect(opnormInf([[1, -2], [3, 0.5]])).toBe(3.5);
    expect(matTvec([[1, 2], [3, 4]], [1, 1])).toEqual([4, 6]);
  });
});
