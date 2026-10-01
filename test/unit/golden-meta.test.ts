// The golden data, the Julia reference project and the port must name the same snapshot
// and the same numerical dependencies (plan §5.1).

import { describe, expect, it } from "vitest";
import { readFileSync } from "node:fs";
import { SNAPSHOT_COMMIT } from "../../src/solver/index.ts";
import { GOLDEN_FILES, decodeNonfinite, loadGolden } from "../golden.ts";

function manifestVersions(): Map<string, string> {
  const text = readFileSync(new URL("../../reference/Manifest.toml", import.meta.url), "utf8");
  const out = new Map<string, string>();
  for (const block of text.split(/^\[\[deps\./m).slice(1)) {
    const name = block.slice(0, block.indexOf("]]"));
    const m = /^version = "([^"]+)"/m.exec(block);
    if (m?.[1]) out.set(name, m[1]);
  }
  return out;
}

describe("golden data provenance", () => {
  const versions = manifestVersions();

  it.each(GOLDEN_FILES)("%s matches the snapshot commit and the Manifest", (name) => {
    const { meta } = loadGolden(name);
    expect(meta.commit).toBe(SNAPSHOT_COMMIT);
    expect(meta.version).toBe(versions.get("ExactMHDRiemannSolver"));
    for (const [pkg, v] of Object.entries(meta.packages)) {
      expect(versions.get(pkg), pkg).toBe(v);
    }
    expect(Object.keys(meta.packages)).toContain("NonlinearSolve");
  });

  it("decodes non-finite numbers", () => {
    const v = JSON.parse('{"a":"Inf","b":"-Inf","c":"NaN","d":"Infinity","e":1.0e-5}', decodeNonfinite);
    expect(v).toEqual({ a: Infinity, b: -Infinity, c: NaN, d: "Infinity", e: 1e-5 });
  });

  it("stores Julia's residual Inf of a refused problem as Infinity", () => {
    const { cases } = loadGolden<{ cases: { name: string; residual: number }[] }>("solve_examples");
    expect(cases.find((c) => c.name === "refusal/gamma")?.residual).toBe(Infinity);
  });
});
