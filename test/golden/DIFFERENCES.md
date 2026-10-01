# Differences between the TS port and the golden data

Every case where the port's `retcode`, `reason` or `method` differs from the golden data,
or where a quantity misses its tolerance of plan §5.3, is listed here with its
explanation (plan §5.3).

| set | case | golden (Julia `6da1c02`) | port | explanation |
| --- | --- | --- | --- | --- |
| `solve_examples` | `refusal/bt_larger_at_threshold` | CheckFailed | Success (expected) | Julia checker false positive on a vanishingly weak fan; the port includes the fix of Julia commit `0ca567a` (`reference/BASELINE.md`). Disappears when the snapshot is re-pinned to the fix. |
| `solve_examples` | `refusal/bt_ok_2` | CheckFailed | Success (expected) | as above |
| `solve_small_bt` | #372, #402, #419, #432, #449 | Success (homotopy) | NoConvergence | a > c_A, larger side Bt/√p ≤ 3e-4: the continuation starts at the degenerate point R = L (singular twist direction); Julia's path winds the twist angle through ~30 turns. FD vs ForwardDiff Jacobians cannot be expected to agree there (reference/BASELINE.md). |
| `solve_small_bt` | #362 | NoConvergence | Success (homotopy) | as above, the other way round |

## Equivalences used in the comparison (`test/compare.ts`)

- **Vanishing waves:** a wave with path variable ψ ≈ ±1e-17 is labelled shock or fan by the sign of round-off (e.g. `single_wave/*`, `refusal/extreme_ratio_at_limit`); same family and zero width count as the same wave.
- **Round-off ties of the canonical frame:** equal Bt/√p on both sides (mirror decision) or an exactly coplanar problem in a rotated frame (twist ±π) can lead to the side-swapped or z-mirrored representation of the same solution; wave tables are compared in the user frame, coplanar ones on mirror-invariant quantities.
- **Ψ** is a diagnostic only (it is not unique under these ties).
