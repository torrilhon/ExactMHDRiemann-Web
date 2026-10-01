// Plan §4.6: in the solver, every sqrt/log/pow goes through the checked wrappers of
// math.ts (sqrtD, logD, powD), which throw DomainError where Julia does, instead of
// returning NaN. This check rejects bare Math.sqrt / Math.log / Math.pow and `**` with a
// non-literal exponent anywhere under src/solver except math.ts itself.
import { readdirSync, readFileSync, statSync } from "node:fs";
import { join, relative } from "node:path";

const ROOT = new URL("../src/solver/", import.meta.url).pathname;
const EXEMPT = new Set(["math.ts"]);
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
  readFileSync(file, "utf8").split("\n").forEach((line, i) => {
    const code = line.replace(/\/\/.*$/, "");
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
