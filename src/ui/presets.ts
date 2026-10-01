// Presets: the examples of the Julia package (examples/*.toml).

import type { ProblemInput } from "../solver/api.ts";

export interface Preset { id: string; label: string; note: string; problem: ProblemInput }

export const PRESETS: Preset[] = [
  {
    id: "paper",
    label: "Torrilhon (2002), Sec. 2.1",
    note: "The example of the report, twist angle 1.5: fast fan, rotation, slow fan, contact, slow shock, rotation, fast shock.",
    problem: {
      gamma: 1.6666666666666667,
      L: [3.0, 0.0, 0.0, 0.0, 1.5, 1.0, 0.0, 3.0],
      R: [1.0, 0.0, 0.0, 0.0, 1.5, 0.0707372016677029, 0.9974949866040544, 1.0],
      t: 0.4, x: [-1.0, 1.0], n: 2001,
    },
  },
  {
    id: "briowu",
    label: "Brio & Wu (1988)",
    note: "Coplanar shock tube, γ = 2, regular solution (the compound wave of numerical schemes is replaced by a rotation and a slow shock).",
    problem: {
      gamma: 2.0,
      L: [1.0, 0.0, 0.0, 0.0, 0.75, 1.0, 0.0, 1.0],
      R: [0.125, 0.0, 0.0, 0.0, 0.75, -1.0, 0.0, 0.1],
      t: 0.1, x: [-0.5, 0.5], n: 2001,
    },
  },
];
