// Fuzz harness (plan §5.4): random problems in the domain of random_problems.csv
// (Bn ∈ [0.3, 2.3], |Bt| ∈ [0.2, 2], ρ, p ∈ [0.2, 5], |vx| ≤ 1, |vt| ≤ 0.5, γ ∈ {1.4, 5/3, 2}),
// plus a share of coplanar and of small-Bt problems. Fails on CheckFailed, on any
// uncaught exception and on a Success whose own checks report messages.
//   npm run fuzz -- [count] [seed]

import { solve } from "../src/solver/solve.ts";
import { riemannProblem } from "../src/solver/types.ts";

const N = Number(process.argv[2] ?? 1000);
let seed = Number(process.argv[3] ?? Date.now() % 2 ** 31);
const seed0 = seed;
function rand(): number {                                    // mulberry32
  seed |= 0; seed = (seed + 0x6d2b79f5) | 0;
  let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
  t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
  return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
}
const U = (a: number, b: number) => a + (b - a) * rand();

function side(Bn: number, smallBt: boolean): number[] {
  const bt = smallBt ? 10 ** U(-4, -2) : U(0.2, 2), th = U(0, 2 * Math.PI);
  const p = U(0.2, 5);
  const b = smallBt ? bt * Math.sqrt(p) : bt;
  return [U(0.2, 5), U(-1, 1), U(-0.5, 0.5), U(-0.5, 0.5), Bn, b * Math.cos(th), b * Math.sin(th), p];
}

const counts = new Map<string, number>();
const bad: string[] = [];
const t0 = performance.now();
for (let i = 0; i < N; i++) {
  const Bn = U(0.3, 2.3);
  const L = side(Bn, false), R = side(Bn, rand() < 0.1);
  if (rand() < 0.2) {                                        // coplanar, antiparallel
    const f = -Math.hypot(R[5]!, R[6]!) / Math.hypot(L[5]!, L[6]!);
    R[5] = f * L[5]!; R[6] = f * L[6]!;
  }
  const gamma = [1.4, 5 / 3, 2][Math.floor(3 * rand())]!;
  const prob = riemannProblem(L, R, gamma);
  try {
    const sol = solve(prob);
    const k = sol.retcode + (sol.reason === "none" ? "" : ":" + sol.reason);
    counts.set(k, (counts.get(k) ?? 0) + 1);
    if (sol.retcode === "CheckFailed" || (sol.retcode === "Success" && (sol.check?.messages.length ?? 0) > 0)) {
      bad.push(`#${i} ${k}: ${sol.check?.messages.join("; ")}\n    ${JSON.stringify({ L, R, gamma })}`);
    }
  } catch (e) {
    bad.push(`#${i} exception ${e instanceof Error ? `${e.name}: ${e.message}` : String(e)}\n    ${JSON.stringify({ L, R, gamma })}`);
  }
}
console.log(`fuzz: ${N} problems, seed ${seed0}, ${((performance.now() - t0) / 1000).toFixed(1)} s`);
console.log("  " + [...counts.entries()].sort((a, b) => b[1] - a[1]).map(([k, v]) => `${k} ${v}`).join(", "));
if (bad.length > 0) {
  console.log(`  FAILURES (${bad.length}):\n  ` + bad.join("\n  "));
  process.exit(1);
}
console.log("  no failures");
