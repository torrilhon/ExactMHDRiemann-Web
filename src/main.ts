import { PORT_VERSION, SNAPSHOT_COMMIT, SNAPSHOT_REPO, SNAPSHOT_VERSION } from "./solver/index.ts";
import "./ui/style.css";

const el = document.querySelector<HTMLElement>("#snapshot");
if (el) {
  const link = document.createElement("a");
  link.href = `${SNAPSHOT_REPO}/tree/${SNAPSHOT_COMMIT}`;
  link.textContent = `ExactMHDRiemannSolver.jl ${SNAPSHOT_VERSION} (${SNAPSHOT_COMMIT.slice(0, 7)})`;
  el.replaceChildren("Port ", PORT_VERSION, " of ", link);
}
