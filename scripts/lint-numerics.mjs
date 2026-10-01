// Plan §4.6: in the solver, every sqrt/log/pow goes through the checked wrappers of
// math.ts (sqrtD, logD, powD), which throw DomainError where Julia does, instead of
// returning NaN. This check rejects bare Math.sqrt / Math.log / Math.pow and `**` with a
// non-literal exponent anywhere under src/solver except math.ts itself (and format.ts, which only formats).
import { readdirSync, readFileSync, statSync } from "node:fs";
import { join, relative } from "node:path";

const ROOT = new URL("../src/solver/", import.meta.url).pathname;
// math.ts defines the wrappers; format.ts and exact.ts are formatting and exact
// (BigInt) arithmetic, not floating-point numerics
// rk/engine.ts: step-size control only (norms, powers and logs of nonnegative error
// estimates); the ODE right-hand sides it calls are checked as usual
const EXEMPT = new Set(["math.ts", "format.ts", "exact.ts", "rk/engine.ts"]);
const RULES = [
  [/\bMath\.(sqrt|log|log2|log10|log1p|pow)\b/, "use sqrtD / logD / powD from math.ts"],
  [/\*\*\s*(?![\d(])/, "use x*x or powD; `**` only with a literal integer exponent"],
];

function* files(dir) {
  for (const name of readdirSync(dir)) {
    const p = join(dir, name);
    if (statSync(p).isDirectory()) yield* files(p);
    else if (p.endsWith(".ts")) yield p;
  }
}

let bad = 0;
for (const file of files(ROOT)) {
  const rel = relative(ROOT, file);
  if (EXEMPT.has(rel)) continue;
  // blank out comments (keeping line numbers), then check the code
  const text = readFileSync(file, "utf8")
    .replace(/\/\*[\s\S]*?\*\//g, (m) => m.replace(/[^\n]/g, " "))
    .replace(/\/\/.*$/gm, "");
  const orig = readFileSync(file, "utf8").split("\n");
  text.split("\n").forEach((code, i) => {
    const line = orig[i] ?? "";
    for (const [re, msg] of RULES) {
      if (re.test(code)) {
        console.error(`src/solver/${rel}:${i + 1}: ${msg}\n    ${line.trim()}`);
        bad++;
      }
    }
  });
}
if (bad > 0) {
  console.error(`\n${bad} unchecked numeric call(s)`);
  process.exit(1);
}
console.log("lint:numerics ok");
