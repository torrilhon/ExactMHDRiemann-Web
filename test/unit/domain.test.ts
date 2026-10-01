import { describe, expect, it } from "vitest";
import { validateInput } from "../../src/solver/canonical.ts";
import { defaultOptions, riemannProblem } from "../../src/solver/types.ts";
import { refusalDetail } from "../../src/ui/domain.ts";

const reason = (L: number[], R: number[]) => validateInput(riemannProblem(L, R, 5 / 3), defaultOptions())[1];

// ρ, vx, vy, vz, Bx, By, Bz, p
const L = [1, 0, 0, 0, 1, 1e-5, 0, 1];
const R = [0.5, 0, 0, 0, 1, 0, 0, 0.25];

describe("refusal details", () => {
  it("names the transverse fields that miss their thresholds", () => {
    expect(reason(L, R)).toBe("switch_on_off");
    expect(refusalDetail("switch_on_off", L, R)).toBe(
      "Measured: larger side (left) |Bt|/√p = 1e-5 < 1e-4; smaller side (right) |Bt|/√p = 0 < 1e-6. " +
      "The regular solver applies because |Bx|/√p = 2 > 1e-10.");
    const R2 = [0.5, 0, 0, 0, 1, 1e-3, 0, 0.25];          // right 2e-3 is the larger side, left 1e-5 suffices
    expect(reason(L, R2)).toBe("none");
  });

  it("gives the jump ratios", () => {
    const R3 = [1e7, 0, 0, 0, 1, 1, 0, 1];
    expect(reason(L, R3)).toBe("extreme_ratio");
    expect(refusalDetail("extreme_ratio", L, R3)).toBe("Measured: ρR/ρL = 1e7, pR/pL = 1; both must lie in [1e-6, 1e6].");
  });

  it("says nothing for other reasons", () => {
    expect(refusalDetail("vacuum", L, R)).toBe("");
  });
});
