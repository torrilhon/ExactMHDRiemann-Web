// Browser smoke test of the built page (plan §5.4): the solver runs in a real browser
// engine (its Math.* may differ from Node's) and the page works end to end.
//   npm run build && npx vite preview --port 4174 & node test/e2e/smoke.mjs [url]

import { chromium } from "playwright";

const url = process.argv[2] ?? "http://localhost:4174/";
const browser = await chromium.launch(process.env.CHROMIUM_PATH ? { executablePath: process.env.CHROMIUM_PATH } : {});
const page = await browser.newPage({ viewport: { width: 1280, height: 900 } });
const errors = [];
page.on("pageerror", (e) => errors.push(String(e)));
page.on("console", (m) => { if (m.type() === "error") errors.push(m.text()); });
await page.goto(url);
const expect = (cond, msg) => { if (!cond) { console.error("FAIL:", msg); process.exitCode = 1; } };

// 1. the report example (computed on load) reproduces Tables 1-2
await page.waitForSelector(".status.good", { timeout: 60_000 });
const cells = await page.$$eval("#waves-body tr", (rows) => rows.map((r) => [...r.children].map((c) => c.textContent)));
const kinds = cells.slice(1).map((r) => r[0]);
expect(JSON.stringify(kinds) === JSON.stringify(["fast fan", "rotation", "slow fan", "contact", "slow shock", "rotation", "fast shock"]), `wave kinds ${kinds}`);
const ref = [-1.474921, -0.631585, -0.521394, 0.402052, 1.279597, 1.568066, 2.072332];
cells.slice(1).forEach((r, i) => expect(Math.abs(Number(r[1]) - ref[i]) <= 2e-6, `s_left of wave ${i + 1}: ${r[1]}`));
expect((await page.locator(".plot").count()) === 8, "8 plots");

// 2. every preset solves
for (const id of ["briowu", "sod"]) {
  await page.selectOption("#preset", id);
  await page.click("#compute");
  await page.waitForFunction(() => !document.querySelector("#compute").disabled);
  await page.waitForTimeout(200);
  expect((await page.getAttribute("#status", "class")).includes("good"), `${id}: ${await page.textContent("#status")}`);
}
expect(errors.length === 0, `console errors: ${errors.join("; ")}`);
console.log(process.exitCode ? "browser smoke test FAILED" : "browser smoke test passed");
await browser.close();
