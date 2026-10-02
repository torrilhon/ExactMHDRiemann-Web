import { describe, expect, it } from "vitest";
import { btAngle, fix, waveTable } from "../../src/ui/wavetable.ts";

describe("wave table", () => {
  it("formats zero without a sign and angles in (-180°, 180°]", () => {
    expect(fix(-1e-17)).toBe("0.000000");
    expect(fix(-0.25)).toBe("-0.250000");
    expect(btAngle(0, 1)).toBe("90.00°");
    expect(btAngle(-1, -0)).toBe("180.00°");
    expect(btAngle(0, 0)).toBe("–");
  });

  it("lists each state once with the waves between them and marks unchanged values", () => {
    const L = [1, 0, 0, 0, 1, 1, 0, 1];
    const lines = waveTable(L, [
      { kind: "contact", s_left: 0.5, s_right: 0.5, rho: 2, vx: 0, vy: 0, vz: 0, Bx: 1, By: 1, Bz: 0, p: 1 },
      { kind: "fast_fan", s_left: 1, s_right: 1.25, rho: 1.5, vx: 0.1, vy: 0, vz: 0, Bx: 1, By: 0.5, Bz: 0, p: 0.5 },
    ]);
    expect(lines.map((l) => l.label)).toEqual(["left state", "contact", "state 1", "fast fan", "right state"]);
    expect(lines[1]).toEqual({ kind: "wave", label: "contact", speed: "s = 0.500000" });
    expect(lines[3]).toEqual({ kind: "wave", label: "fast fan", speed: "s = 1.000000 … 1.250000" });
    const s1 = lines[2]!;
    expect(s1.kind === "state" && s1.same).toEqual([false, true, true, true, true, true, true, true]);
    expect(s1.kind === "state" && s1.cells[6]).toBe("0.00°");
  });
});
