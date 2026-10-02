// Small multiples of the 8 primitive variables over x (uPlot), one series each, with a
// crosshair synchronized across the panels and the value under the cursor in each
// panel's header. Colors come from CSS custom properties, re-read on theme changes.

import uPlot from "uplot";
import "uplot/dist/uPlot.min.css";

let charts: uPlot[] = [];
let lastData: { x: number[]; W: number[][]; names: string[]; host: HTMLElement } | null = null;

const css = (name: string) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

function panelWidth(host: HTMLElement): number {
  const cols = host.clientWidth >= 720 ? 2 : 1;
  const gap = 16;
  return Math.max(240, Math.floor((host.clientWidth - gap * (cols - 1)) / cols));
}

export function renderPlots(host: HTMLElement, x: number[], W: number[][], names: string[]): void {
  lastData = { x, W, names, host };
  for (const c of charts) c.destroy();
  charts = [];
  host.replaceChildren();
  const ink = css("--muted"), grid = css("--grid"), line = css("--series-1");
  const sync = uPlot.sync("profiles");
  const width = panelWidth(host);
  names.forEach((name, j) => {
    const fig = document.createElement("figure");
    fig.className = "plot";
    const cap = document.createElement("figcaption");
    const title = document.createElement("span");
    title.className = "plot-title";
    title.textContent = name;
    const val = document.createElement("span");
    val.className = "plot-value";
    cap.append(title, val);
    fig.append(cap);
    host.append(fig);
    const y = W.map((row) => row[j]!);
    const opts: uPlot.Options = {
      width, height: 170,
      legend: { show: false },
      cursor: { sync: { key: sync.key }, points: { size: 8 }, drag: { x: false, y: false } },
      scales: { x: { time: false } },
      axes: [
        { stroke: ink, grid: { stroke: grid, width: 1 }, ticks: { stroke: grid, width: 1 }, size: 32 },
        { stroke: ink, grid: { stroke: grid, width: 1 }, ticks: { stroke: grid, width: 1 }, size: 56 },
      ],
      // pxAlign 0: no snapping to device pixels, which on fractional pixel ratios can split
      // values equal up to round-off (1e-16) into neighbouring rows, a step that is not there
      series: [{}, { label: name, stroke: line, width: 2, points: { show: false }, pxAlign: 0 }],
      hooks: {
        setCursor: [(u) => {
          const i = u.cursor.idx;
          val.textContent = i == null ? "" : `x = ${x[i]!.toFixed(4)}   ${name} = ${y[i]!.toPrecision(6)}`;
        }],
      },
    };
    charts.push(new uPlot(opts, [x, y], fig));
  });
}

function rerender(): void {
  if (lastData) renderPlots(lastData.host, lastData.x, lastData.W, lastData.names);
}

let resizeTimer: ReturnType<typeof setTimeout> | undefined;
addEventListener("resize", () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(() => {
    if (!lastData) return;
    const w = panelWidth(lastData.host);
    for (const c of charts) c.setSize({ width: w, height: 170 });
  }, 100);
});
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", rerender);
