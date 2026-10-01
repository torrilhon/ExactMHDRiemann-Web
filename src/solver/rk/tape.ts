// Step tape for frozen-step finite differences (plan §4.3). For the base point of a
// Jacobian the tape records the step sequence of each fan integration, keyed by its slot
// (fast/slow fan of the left/right side; a residual has at most one of each); for the
// perturbed points it replays them, so that the perturbed residuals come from the same
// discrete map (internal numerical differentiation) and the difference quotients are
// free of step-size-control noise.

import type { OdeOptions, OdeSolution, Rhs, Tableau } from "./engine.ts";
import { integrate } from "./engine.ts";

export class StepTape {
  mode: "off" | "record" | "replay" = "off";
  private seqs = new Map<string, number[]>();

  record(): void {
    this.mode = "record";
    this.seqs.clear();
  }

  replay(): void {
    this.mode = "replay";
  }

  off(): void {
    this.mode = "off";
  }

  /** The step sequence recorded for a slot. */
  recorded(slot: string): number[] | undefined {
    return this.seqs.get(slot);
  }

  /** Store a step sequence for a slot (record mode), e.g. from a cached integration. */
  put(slot: string, steps: number[]): void {
    if (this.mode === "record") this.seqs.set(slot, steps);
  }

  /** `integrate` through the tape. In replay mode a missing sequence falls back to adaptive. */
  integrate(slot: string, tab: Tableau, f: Rhs, y0: readonly number[], t0: number, t1: number,
    opts: OdeOptions): OdeSolution {
    if (this.mode === "replay") {
      const steps = this.seqs.get(slot);
      if (steps !== undefined) return integrate(tab, f, y0, t0, t1, { ...opts, steps });
      return integrate(tab, f, y0, t0, t1, opts);
    }
    const sol = integrate(tab, f, y0, t0, t1, opts);
    if (this.mode === "record") this.seqs.set(slot, sol.steps);
    return sol;
  }
}

/**
 * Memo of fan end states for one solve: the residual's fans depend only on their outer
 * state and path variable, and a finite-difference Jacobian perturbs one component at a
 * time, so most of its fan integrations repeat earlier ones exactly. Entries keep their
 * step sequence so that a cache hit can still feed the step tape. Small LRU.
 */
export class FanCache {
  private map = new Map<string, { end: number[]; steps: number[] }>();
  readonly capacity: number;
  constructor(capacity = 32) {
    this.capacity = capacity;
  }

  get(key: string): { end: number[]; steps: number[] } | undefined {
    const v = this.map.get(key);
    if (v !== undefined) { this.map.delete(key); this.map.set(key, v); }
    return v;
  }

  set(key: string, v: { end: number[]; steps: number[] }): void {
    this.map.set(key, v);
    if (this.map.size > this.capacity) this.map.delete(this.map.keys().next().value!);
  }
}
