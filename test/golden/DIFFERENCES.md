# Differences between the TS port and the golden data

Every case where the port's `retcode`, `reason` or `method` differs from the golden data,
or where a quantity misses its tolerance of plan §5.3, is listed here with its
explanation (plan §5.3).

| set | case | golden (Julia `6da1c02`) | port | explanation |
| --- | --- | --- | --- | --- |
| `solve_examples` | `refusal/bt_larger_at_threshold` | CheckFailed | Success (expected) | Julia checker false positive on a vanishingly weak fan; the port includes the fix of Julia commit `0ca567a` (`reference/BASELINE.md`). Disappears when the snapshot is re-pinned to the fix. |
| `solve_examples` | `refusal/bt_ok_2` | CheckFailed | Success (expected) | as above |
