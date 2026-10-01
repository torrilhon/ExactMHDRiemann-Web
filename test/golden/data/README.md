# Test data

- `random_problems.csv`: 2000 random non-planar problems (γ = 5/3). Columns: seed, L (8), R (8),
  primitive states (ρ, vx, vy, vz, Bx, By, Bz, p). Generated with the seeded generator of the
  2026 code review (Bn ∈ [0.3, 2.3], |Bt| ∈ [0.2, 2], ρ, p ∈ [0.2, 5], |vx| ≤ 1, |vt| ≤ 0.5).
- `c_reference.csv`: for 355 of the first 400 problems, the solution of the original C code
  (with the root-finder and NaN fixes of the review) that passed the independent Python checks.
  Columns: seed, L (8), R (8), state left of the contact (8), state right of the contact (8).
  Values carry 8 decimals, the C code's output precision.
