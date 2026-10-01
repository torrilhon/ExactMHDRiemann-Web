import { describe, expect, it } from "vitest";
import { readFileSync } from "node:fs";
import { SNAPSHOT_COMMIT } from "../../src/solver/index.ts";

describe("snapshot provenance", () => {
  it("is a full commit hash", () => {
    expect(SNAPSHOT_COMMIT).toMatch(/^[0-9a-f]{40}$/);
  });

  it("matches the Julia reference project", () => {
    // reference/Manifest.toml pins ExactMHDRiemannSolver to the snapshot commit
    const manifest = readFileSync(new URL("../../reference/Manifest.toml", import.meta.url), "utf8");
    const block = manifest.split("[[deps.ExactMHDRiemannSolver]]")[1]?.split("[[deps.")[0] ?? "";
    expect(block).toContain(`repo-rev = "${SNAPSHOT_COMMIT}"`);
  });
});
