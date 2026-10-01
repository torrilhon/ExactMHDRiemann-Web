import { describe, expect, it } from "vitest";
import { readFileSync } from "node:fs";
import { problemLoad, problemToml } from "../../src/solver/problemfile.ts";
import { ArgumentError } from "../../src/solver/math.ts";

// examples/ of the Julia package, copied by reference/generate_golden.jl
const ex = (f: string) => readFileSync(new URL(`../golden/examples/${f}`, import.meta.url), "utf8");

describe("problem files (problem_load)", () => {
  it("reads the examples of the Julia package as Julia does", () => {
    const p = problemLoad(ex("paper.toml"), "examples/paper.toml");
    expect(p.L).toEqual([3, 0, 0, 0, 1.5, 1, 0, 3]);
    expect(p.gamma).toBeCloseTo(5 / 3, 15);
    expect([p.t, ...p.x, p.n]).toEqual([0.4, -1, 1, 2001]);
    const b = problemLoad(ex("briowu.toml"), "briowu.toml");
    expect(b.gamma).toBe(2);
    expect(b.R).toEqual([0.125, 0, 0, 0, 0.75, -1, 0, 0.1]);
    expect(b.x[0]).toBe(-0.5);
  });
  it("applies Julia's defaults", () => {
    const p = problemLoad("left = [1, 0, 0, 0, 1, 1, 0, 1]\nright = [1, 0, 0, 0, 1, 1, 0, 1]\n", "f.toml");
    expect(p.gamma).toBe(5 / 3);
    expect([p.t, ...p.x, p.n]).toEqual([1, -1, 1, 2001]);
  });
  it("gives Julia's error messages", () => {
    const err = (s: string) => { try { problemLoad(s, "f.toml"); return ""; } catch (e) { expect(e).toBeInstanceOf(ArgumentError); return (e as Error).message; } };
    expect(err("left = [1, 0, 0]\nright = [1, 0, 0, 0, 1, 1, 0, 1]\n")).toBe("f.toml: `left` must be 8 numbers (ρ, vx, vy, vz, Bx, By, Bz, p)");
    expect(err("left = [1, 0, 0, 0, 1, 1, 0, 1]\n")).toBe("f.toml: missing `right`");
    expect(err("left = [1, 0, 0, 0, 1, 1, 0, 1]\nright = [1, 0, 0, 0, 1, 1, 0, 1]\n[output]\nx = [1]\n")).toBe("f.toml: `output.x` must be [xmin, xmax]");
    expect(err("left = [1, ")).toMatch(/^f\.toml: not a valid TOML file/);
  });
  it("round-trips through problemToml", () => {
    const p = { L: [3, 0, 0, 0, 1.5, 1, 0, 3], R: [1, 0, 0, 0, 1.5, 0.0707372016677029, 0.9974949866040544, 1], gamma: 5 / 3, t: 0.4, x: [-1, 1] as [number, number], n: 2001 };
    const text = problemToml(p);
    expect(problemLoad(text, "problem.toml")).toEqual(p);
    expect(text).toContain("gamma = 1.6666666666666667");
    expect(text).toContain("left  = [3.0, 0.0, 0.0, 0.0, 1.5, 1.0, 0.0, 3.0]");
  });
});
