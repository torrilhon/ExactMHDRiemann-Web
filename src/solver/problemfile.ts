// Problem files (TOML), the format of Julia's `problem_load` (problemfile.jl):
//
//   gamma = 1.6666666666666667            # optional, default 5/3
//   left  = [rho, vx, vy, vz, Bx, By, Bz, p]
//   right = [rho, vx, vy, vz, Bx, By, Bz, p]
//   [output]                              # optional
//   t = 0.4                               # default 1.0
//   x = [-1.0, 1.0]                       # default [-1, 1]
//   n = 2001                              # default 2001

import { parse } from "smol-toml";
import { juliaRepr } from "./format.ts";
import { ArgumentError } from "./math.ts";
import type { ProblemInput } from "./api.ts";
import { PORT_VERSION, SNAPSHOT_COMMIT } from "./snapshot.ts";

// Julia: `x isa Real` (Bool is a Real there, true → 1.0)
const isReal = (x: unknown): x is number | bigint | boolean =>
  typeof x === "number" || typeof x === "bigint" || typeof x === "boolean";
const toFloat = (x: number | bigint | boolean) => (typeof x === "boolean" ? (x ? 1 : 0) : Number(x));

/**
 * Julia `problem_load` on the text of a file; `path` is used in the error messages
 * (ArgumentError, as in Julia). TOML syntax errors are reported as ArgumentError too.
 */
export function problemLoad(text: string, path: string): ProblemInput {
  let cfg: Record<string, unknown>;
  try {
    cfg = parse(text) as Record<string, unknown>;
  } catch (e) {
    throw new ArgumentError(`${path}: not a valid TOML file (${e instanceof Error ? e.message.split("\n")[0] : String(e)})`);
  }
  const state = (key: string): number[] => {
    if (!(key in cfg)) throw new ArgumentError(`${path}: missing \`${key}\``);
    const v = cfg[key];
    if (!(Array.isArray(v) && v.length === 8 && v.every(isReal))) {
      throw new ArgumentError(`${path}: \`${key}\` must be 8 numbers (ρ, vx, vy, vz, Bx, By, Bz, p)`);
    }
    return v.map(toFloat);
  };
  const L = state("left"), R = state("right");
  const num = (v: unknown, what: string): number => {
    if (!isReal(v)) throw new ArgumentError(`${path}: \`${what}\` must be a number`);
    return toFloat(v);
  };
  const gamma = "gamma" in cfg ? num(cfg.gamma, "gamma") : 5 / 3;
  const out = (cfg.output ?? {}) as Record<string, unknown>;
  const t = "t" in out ? num(out.t, "output.t") : 1.0;
  const xr = "x" in out ? out.x : [-1.0, 1.0];
  if (!(Array.isArray(xr) && xr.length === 2 && xr.every(isReal))) {
    throw new ArgumentError(`${path}: \`output.x\` must be [xmin, xmax]`);
  }
  const n = "n" in out ? num(out.n, "output.n") : 2001;
  if (!Number.isInteger(n)) throw new ArgumentError(`${path}: \`output.n\` must be an integer`);    // Julia: InexactError
  const [x0, x1] = xr as [number | bigint | boolean, number | bigint | boolean];
  return { L, R, gamma, t, x: [toFloat(x0), toFloat(x1)], n };
}

const vec = (v: readonly number[]) => "[" + v.map(juliaRepr).join(", ") + "]";

/** The input as a problem file, readable by Julia's `problem_load` and `bin/riemann.jl`. */
export function problemToml(p: ProblemInput): string {
  return [
    `# Riemann problem exported by ExactMHDRiemann-Web ${PORT_VERSION} (snapshot ${SNAPSHOT_COMMIT.slice(0, 7)}).`,
    "# Reproduce in Julia:  L, R, γ, t, x = problem_load(\"problem.toml\"); sol = solve(RiemannProblem(L, R; γ))",
    `gamma = ${juliaRepr(p.gamma)}`,
    `left  = ${vec(p.L)}`,
    `right = ${vec(p.R)}`,
    "",
    "[output]",
    `t = ${juliaRepr(p.t)}`,
    `x = ${vec(p.x)}`,
    `n = ${p.n}`,
    "",
  ].join("\n");
}
