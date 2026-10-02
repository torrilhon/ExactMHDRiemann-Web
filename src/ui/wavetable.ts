// The wave table of the page: each state once, from left to right, with a row for the
// wave between each pair of neighbouring states. Pure formatting, no DOM.

import type { WaveRow } from "../solver/output.ts";

export const STATE_COLUMNS = ["ρ", "vx", "vy", "vz", "By", "Bz", "∠Bt", "p"] as const;

export interface StateLine { kind: "state"; label: string; cells: string[]; same: boolean[] }
export interface WaveLine { kind: "wave"; label: string; speed: string }
export type TableLine = StateLine | WaveLine;

/** Fixed six decimals; a value that rounds to zero is shown without a sign. */
export const fix = (v: number, d = 6): string => {
  if (!Number.isFinite(v)) return String(v);
  const s = v.toFixed(d);
  return Number(s) === 0 ? (0).toFixed(d) : s;
};

/** Direction of Bt in degrees, atan2(Bz, By) in (-180°, 180°]; "–" for Bt = 0. */
export function btAngle(By: number, Bz: number): string {
  if (Math.hypot(By, Bz) <= 1e-12) return "–";
  const a = (Math.atan2(Bz, By) * 180) / Math.PI;
  const s = fix(a, 2);
  return `${s === "-180.00" ? "180.00" : s}°`;
}

const stateCells = (W: readonly number[]): string[] =>
  [fix(W[0]!), fix(W[1]!), fix(W[2]!), fix(W[3]!), fix(W[5]!), fix(W[6]!), btAngle(W[5]!, W[6]!), fix(W[7]!)];

const speed = (w: WaveRow): string =>
  fix(w.s_left) === fix(w.s_right) ? `s = ${fix(w.s_left)}` : `s = ${fix(w.s_left)} … ${fix(w.s_right)}`;

/** Lines of the table for the left state L and the waves (each with the state to its right). */
export function waveTable(L: readonly number[], waves: readonly WaveRow[]): TableLine[] {
  const states = [L, ...waves.map((w) => [w.rho, w.vx, w.vy, w.vz, w.Bx, w.By, w.Bz, w.p])];
  const lines: TableLine[] = [];
  let prev: string[] | null = null;
  states.forEach((W, k) => {
    const cells = stateCells(W);
    const label = k === 0 ? "left state" : k === states.length - 1 ? "right state" : `state ${k}`;
    lines.push({ kind: "state", label, cells, same: cells.map((c, i) => prev !== null && prev[i] === c) });
    prev = cells;
    const w = waves[k];
    if (w) lines.push({ kind: "wave", label: w.kind.replace("_", " "), speed: speed(w) });
  });
  return lines;
}
