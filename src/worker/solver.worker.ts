// Web Worker running the solver off the main thread.
// Message in: { id, problem: ProblemInput } → out: { id, ok: true, result } | { id, ok: false, error }.

import { runProblem, type ProblemInput } from "../solver/api.ts";

self.onmessage = (ev: MessageEvent<{ id: number; problem: ProblemInput }>) => {
  const { id, problem } = ev.data;
  try {
    self.postMessage({ id, ok: true, result: runProblem(problem) });
  } catch (e) {
    self.postMessage({ id, ok: false, error: e instanceof Error ? `${e.name}: ${e.message}` : String(e) });
  }
};
