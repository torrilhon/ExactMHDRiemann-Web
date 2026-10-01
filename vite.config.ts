import { defineConfig } from "vitest/config";

export default defineConfig({
  // relative asset paths, so the build works under https://<owner>.github.io/<repo>/
  base: "./",
  worker: { format: "es" },
  test: {
    include: ["test/**/*.test.ts"],
  },
});
