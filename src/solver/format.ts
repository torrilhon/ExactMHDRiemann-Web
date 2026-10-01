// Julia-compatible output formatting (plan §4.7): `repr(::Float64)` and the points of
// `range(a, b; length = n)`, so that the CSV export agrees with Julia's `write_csv`
// character for character.

import { exactRational, ratToDouble } from "./exact.ts";

/**
 * Julia `repr(x)` for a Float64: the shortest round-trip digits (as Ryu), in fixed
 * notation when the decimal point position pt satisfies -4 < pt <= 6, otherwise
 * `d.ddde±n`; always with a decimal point; `-0.0`, `NaN`, `Inf`, `-Inf`.
 * (Julia base/ryu/shortest.jl, `writeshortest` with precision = -1.)
 */
export function juliaRepr(x: number): string {
  if (Number.isNaN(x)) return "NaN";
  if (x === Infinity) return "Inf";
  if (x === -Infinity) return "-Inf";
  if (x === 0) return Object.is(x, -0) ? "-0.0" : "0.0";
  const neg = x < 0;
  // shortest round-trip digits: "d.ddde±n"
  const [mant = "", exps = "0"] = Math.abs(x).toExponential().split("e");
  const digits = mant.replace(".", "");
  const olength = digits.length;
  const e = Number(exps);
  const pt = e + 1;                          // value = 0.d1d2… × 10^pt
  const nexp = pt - olength;                 // value = digits × 10^nexp
  let expForm = true;
  if (-4 < pt && pt <= 6) {
    // Julia also requires !(pt >= olength && |mod(x + 0.05, 10^(pt - olength)) - 0.05| > 0.05)
    const bad = pt >= olength && Math.abs(juliaMod(x + 0.05, 10 ** nexp) - 0.05) > 0.05;
    if (!bad) expForm = false;
  }
  let s: string;
  if (!expForm) {
    if (pt <= 0) s = "0." + "0".repeat(-pt) + digits;
    else if (pt >= olength) s = digits + "0".repeat(pt - olength) + ".0";
    else s = digits.slice(0, pt) + "." + digits.slice(pt);
  } else {
    const rest = digits.slice(1);
    s = (digits[0] ?? "") + "." + (rest === "" ? "0" : rest) + "e" + String(pt - 1);
  }
  return neg ? "-" + s : s;
}

/** Julia `mod(a, b)` for floats: result has the sign of b. */
function juliaMod(a: number, b: number): number {
  const r = a % b;
  return r !== 0 && (r < 0) !== (b < 0) ? r + b : r;
}

// ---------------------------------------------------------------- ranges

const MAXINT_F32 = 16777216;               // maxintfloat(Float32, Int)
const MAXINT_F64 = 9007199254740992;       // maxintfloat(Float64, Int)

/** Julia `rat(x)` (base/twiceprecision.jl): best rational approximation with |a|, |b| ≤ 2^24. */
export function juliaRat(x: number): [number, number] {
  let y = x;
  let a = 1, d = 1, b = 0, c = 0;
  while (Math.abs(y) <= MAXINT_F32) {
    const f = Math.trunc(y);
    y -= f;
    [a, c] = [f * a + c, a];
    [b, d] = [f * b + d, b];
    if (Math.max(Math.abs(a), Math.abs(b)) > MAXINT_F32) return [c, d];
    if (a / b === x) break;
    y = 1 / y;
  }
  return [a, b];
}

function gcd(a: number, b: number): number {
  a = Math.abs(a); b = Math.abs(b);
  while (b !== 0) [a, b] = [b, a % b];
  return a;
}

/**
 * Points of Julia's `range(start, stop; length = len)` for Float64.
 *
 * Julia represents the endpoints as rationals start_n/den, stop_n/den (via `rat`) when
 * that is exact in Float64 and evaluates ((len - i) start_n + (i - 1) stop_n) / ((len - 1) den)
 * in double-double precision; this is reproduced with exact BigInt arithmetic and correct
 * rounding. Otherwise the exact binary values of the endpoints are interpolated (Julia's
 * double-double path agrees with that up to rare last-bit cases).
 */
export function juliaRange(start: number, stop: number, len: number): number[] {
  if (!Number.isInteger(len) || len < 0) throw new RangeError(`range: invalid length ${len}`);
  if (len < 2) {
    if (len === 1 && start !== stop) {
      throw new RangeError(`range(${start}, stop=${stop}, length=${len}): endpoints differ`);
    }
    return len === 1 ? [start] : [];
  }
  if (start === stop) return new Array<number>(len).fill(start);
  const [sn0, sd] = juliaRat(start);
  const [en0, ed] = juliaRat(stop);
  if (sd !== 0 && ed !== 0) {
    const den = sd * (ed / gcd(sd, ed));              // lcm_unchecked
    if (den !== 0 && Math.abs(den * start) <= MAXINT_F64 && Math.abs(den * stop) <= MAXINT_F64) {
      const sn = roundEven(den * start);              // round(Int, ·): RoundNearest, ties to even
      const en = roundEven(den * stop);
      if (sn / den === start && en / den === stop) {
        return rationalPoints(BigInt(sn), BigInt(en), BigInt(den), len);
      }
    }
  }
  const [an, ad] = exactRational(start);
  const [bn, bd] = exactRational(stop);
  // start = an/ad, stop = bn/bd → common denominator
  return rationalPoints(an * bd, bn * ad, ad * bd, len);
}

function roundEven(x: number): number {
  const r = Math.round(x);                            // ties towards +Inf
  return r - x === 0.5 && r % 2 !== 0 ? r - 1 : r;
}

function rationalPoints(sn: bigint, en: bigint, den: bigint, len: number): number[] {
  const out = new Array<number>(len);
  const D = BigInt(len - 1) * den;
  for (let i = 1; i <= len; i++) {
    const N = BigInt(len - i) * sn + BigInt(i - 1) * en;
    out[i - 1] = ratToDouble(N, D);
  }
  return out;
}
