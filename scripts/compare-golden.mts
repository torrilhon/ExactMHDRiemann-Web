// Runs the TS solver on a golden solve set and reports agreement with Julia (plan §5.3).
//   node --experimental-strip-types scripts/compare-golden.mts solve_stress [limit]
// Writes a summary to stdout and the per-case results to test/golden/report/<set>.json.

import { mkdirSync, writeFileSync } from "node:fs";
import { caseLabel, runCase, type GoldenCase } from "../test/compare.ts";
import { loadGolden, type GoldenName } from "../test/golden.ts";

const set = (process.argv[2] ?? "solve_examples") as GoldenName;
const limit = Number(process.argv[3] ?? Infinity);
const G = loadGolden<{ cases?: GoldenCase[]; scan?: GoldenCase[] }>(set);
const cases = [...(G.cases ?? []), ...(G.scan ?? [])].slice(0, limit);

let codeOk = 0, methodOnly = 0, compared = 0, tableOk = 0, coplanar = 0, checkBad = 0, worstTable = 0, worstPsi = 0;
let ms = 0, julia = 0, tmax = 0;
const mismatches: string[] = [], tableBad: string[] = [], checkFails: string[] = [];
const out: unknown[] = [];
const started = performance.now();
cases.forEach((c, i) => {
  const r = runCase(c, caseLabel(c, i));
  ms += r.ms; julia += c.time_ms; tmax = Math.max(tmax, r.ms);
  if (r.codeMatch) codeOk++;
  else mismatches.push(`${r.label}: julia ${c.retcode}/${c.reason}/${c.method}  ts ${r.sol.retcode}/${r.sol.reason}/${r.sol.method}`);
  // a different method (direct vs homotopy) with the same outcome still compares the solution
  const sameOutcome = r.sol.retcode === c.retcode && r.sol.reason === c.reason;
  if (!r.codeMatch && sameOutcome) methodOnly++;
  if (sameOutcome) {
    compared++;
    if (r.tableErr <= r.tableTol) tableOk++;
    else tableBad.push(`${r.label}: table relerr ${r.tableErr.toExponential(2)} (tol ${r.tableTol.toExponential(1)})`);
    if (Number.isFinite(r.tableErr)) worstTable = Math.max(worstTable, r.tableErr);
    if (r.sol.retcode === "Success" && !r.coplanar) worstPsi = Math.max(worstPsi, r.psiErr);
  }
  if (r.coplanar) coplanar++;
  if (r.sol.retcode === "CheckFailed" || (r.sol.retcode === "Success" && !r.checkOk)) {
    checkBad++;
    checkFails.push(`${r.label}: ${r.messages.join("; ")}`);
  }
  out.push({ label: r.label, julia: [c.retcode, c.reason, c.method], ts: [r.sol.retcode, r.sol.reason, r.sol.method],
    tableErr: r.tableErr, coplanar: r.coplanar, psiErr: r.psiErr, ms: r.ms, julia_ms: c.time_ms,
    maxerr: r.sol.check?.maxerr ?? null });
  if ((i + 1) % 100 === 0) process.stderr.write(`  ${i + 1}/${cases.length} (${((performance.now() - started) / 1000).toFixed(0)} s)\n`);
});
const n = cases.length;
console.log(`${set}: ${n} problems`);
console.log(`  retcode/reason/method identical: ${codeOk}/${n} (${(100 * codeOk / n).toFixed(2)} %)`);
console.log(`  same retcode and reason, different method: ${methodOnly}`);
console.log(`  wave tables within tolerance:    ${tableOk}/${compared}  (worst ${worstTable.toExponential(2)}; coplanar mirror cases ${coplanar})`);
console.log(`  Ψ worst relerr (non-coplanar Success): ${worstPsi.toExponential(2)}`);
console.log(`  CheckFailed or Success with messages: ${checkBad}`);
console.log(`  time: TS total ${(ms / 1000).toFixed(1)} s, median-ish mean ${(ms / n).toFixed(1)} ms, max ${tmax.toFixed(0)} ms; Julia total ${(julia / 1000).toFixed(1)} s`);
for (const [title, list] of [["retcode mismatches", mismatches], ["table mismatches", tableBad], ["check failures", checkFails]] as const) {
  if (list.length > 0) console.log(`  ${title}:\n    ` + list.slice(0, 40).join("\n    ") + (list.length > 40 ? `\n    … ${list.length - 40} more` : ""));
}
mkdirSync(new URL("../test/golden/report/", import.meta.url), { recursive: true });
writeFileSync(new URL(`../test/golden/report/${set}.json`, import.meta.url), JSON.stringify(out) + "\n");
