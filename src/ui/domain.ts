// Measured values behind a refused input: the quantities `validateInput` compares with
// its thresholds, so the page can say by how much a problem misses the supported domain.

import { defaultOptions } from "../solver/types.ts";

const D = defaultOptions();
const num = (v: number) => (v === 0 ? "0" : Math.abs(v) >= 1e-3 && Math.abs(v) < 1e4
  ? String(Number(v.toPrecision(3)))
  : Number(v.toPrecision(2)).toExponential().replace("e+", "e"));
const btp = (W: readonly number[]) => Math.hypot(W[5]!, W[6]!) / Math.sqrt(W[7]!);
const bnp = (W: readonly number[]) => Math.abs(W[4]!) / Math.sqrt(W[7]!);

/** One sentence with the measured values for a refusal reason, or "" for any other reason. */
export function refusalDetail(reason: string, L: readonly number[], R: readonly number[], bnEuler = D.bn_euler): string {
  if (reason === "switch_on_off") {
    const bl = btp(L), br = btp(R);
    const [big, small] = bl >= br ? ["left", "right"] : ["right", "left"];
    const vbig = Math.max(bl, br), vsmall = Math.min(bl, br);
    const miss: string[] = [];
    if (!(vbig >= (1 - 1e-12) * D.bt_min_larger)) miss.push(`larger side (${big}) |Bt|/√p = ${num(vbig)} < ${num(D.bt_min_larger)}`);
    if (!(vsmall >= (1 - 1e-12) * D.bt_min_smaller)) miss.push(`smaller side (${small}) |Bt|/√p = ${num(vsmall)} < ${num(D.bt_min_smaller)}`);
    const bn = Math.max(bnp(L), bnp(R));
    return `Measured: ${miss.join("; ")}. The regular solver applies because |Bx|/√p = ${num(bn)} > ${num(bnEuler)}.`;
  }
  if (reason === "extreme_ratio") {
    const lo = num(1 / D.ratio_max), hi = num(D.ratio_max);
    return `Measured: ρR/ρL = ${num(R[0]! / L[0]!)}, pR/pL = ${num(R[7]! / L[7]!)}; both must lie in [${lo}, ${hi}].`;
  }
  if (reason === "bn_jump") {
    return `Measured: Bx = ${L[4]} on the left, ${R[4]} on the right.`;
  }
  return "";
}
