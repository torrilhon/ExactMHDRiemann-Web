// Step tape for frozen-step finite differences (plan §4.3). The residual evaluates up to
// four fan ODEs. For the base point of a Jacobian the tape records the step sequence of
// every integration in call order; for the perturbed points it replays them, so that the
// perturbed residuals come from the same discrete map (internal numerical
// differentiation) and the difference quotients are free of step-size-control noise.

import type { OdeOptions, OdeSolution, Rhs, Tableau } from "./engine.ts";
import { integrate } from "./engine.ts";

export class StepTape {
  mode: "off" | "record" | "replay" = "off";
  private seqs: number[][] = [];
  private next = 0;

  record(): void {
    this.mode = "record";
    this.seqs = [];
    this.next = 0;
  }

  replay(): void {
    this.mode = "replay";
    this.next = 0;
  }

  off(): void {
    this.mode = "off";
  }

  /** `integrate` through the tape. In replay mode a missing sequence falls back to adaptive. */
  integrate(tab: Tableau, f: Rhs, y0: readonly number[], t0: number, t1: number, opts: OdeOptions): OdeSolution {
    if (this.mode === "replay") {
      const steps = this.seqs[this.next++];
      if (steps !== undefined) return integrate(tab, f, y0, t0, t1, { ...opts, steps });
      return integrate(tab, f, y0, t0, t1, opts);
    }
    const sol = integrate(tab, f, y0, t0, t1, opts);
    if (this.mode === "record") this.seqs.push(sol.steps);
    return sol;
  }
}
