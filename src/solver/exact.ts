// Exact arithmetic helpers (BigInt): exact rational value of a double and correctly
// rounded conversion back. Used by the Julia-compatible range (format.ts) and as the
// fallback of the FMA emulation (math.ts).

/** Exact binary value of a finite double as numerator/denominator (denominator a power of 2). */
export function exactRational(x: number): [bigint, bigint] {
  if (x === 0) return [0n, 1n];
  const buf = new DataView(new ArrayBuffer(8));
  buf.setFloat64(0, x);
  const hi = buf.getUint32(0), lo = buf.getUint32(4);
  const neg = hi >>> 31 === 1;
  const bexp = (hi >>> 20) & 0x7ff;
  let mant = (BigInt(hi & 0xfffff) << 32n) | BigInt(lo);
  let e: number;
  if (bexp === 0) e = -1074;
  else { mant |= 1n << 52n; e = bexp - 1075; }
  if (neg) mant = -mant;
  return e >= 0 ? [mant << BigInt(e), 1n] : [mant, 1n << BigInt(-e)];
}

/** Correctly rounded (round-half-even) Float64 value of the rational N/D. */
export function ratToDouble(N: bigint, D: bigint): number {
  if (D < 0n) { N = -N; D = -D; }
  if (N === 0n) return 0;
  const LIM = 1n << 53n;
  const neg = N < 0n;
  let n = neg ? -N : N;
  if (n <= LIM && D <= LIM) {
    const v = Number(n) / Number(D);        // both exact: IEEE division is correctly rounded
    return neg ? -v : v;
  }
  // scale so that the integer quotient has 54..55 significant bits
  const bitlen = (v: bigint) => v.toString(2).length;
  let shift = 55 - (bitlen(n) - bitlen(D));
  let num = shift >= 0 ? n << BigInt(shift) : n;
  let den = shift >= 0 ? D : D << BigInt(-shift);
  let q = num / den;
  let r = num % den;
  while (bitlen(q) > 55) { den <<= 1n; shift--; q = num / den; r = num % den; }
  while (bitlen(q) < 54) { num <<= 1n; shift++; q = num / den; r = num % den; }
  // value = (q + r/den) × 2^-shift with q of 54 or 55 bits; keep 53 bits, round half even
  let extra = bitlen(q) - 53;
  // subnormal range: keep fewer bits
  const e2 = bitlen(q) - 1 - shift;         // exponent of the leading bit
  if (e2 < -1022) extra += -1022 - e2;
  const mask = (1n << BigInt(extra)) - 1n;
  const rem = q & mask;
  let m = q >> BigInt(extra);
  const half = 1n << BigInt(extra - 1);
  if (rem > half || (rem === half && (r > 0n || (m & 1n) === 1n))) m += 1n;
  const v = Number(m) * Math.pow(2, extra - shift);
  return neg ? -v : v;
}
