# Derived Relations

Formulas derived by combining two or more formulas from
[docs/ANALYTICAL_FORMULAS.md](../../ANALYTICAL_FORMULAS.md).
Each derivation is numerically verified against Liouvillian
eigenvalues (N=2-5).

## Index

| File | Derives | From formulas | Status |
|------|---------|---------------|--------|
| [D01](D01_BANDWIDTH_MODE_DENSITY.md) | BW = 8J cos(pi/N) | 2 | VERIFIED |
| [D02](D02_VEFFECT_QMAX_QMEAN.md) | V(N) = Q_max/Q_mean | 6 + 7 | VERIFIED |
| [D03](D03_CROSSING_TIME_RATIOS.md) | t_X/t_Z = 2.320 | 27 | VERIFIED |
| [D04](D04_DIMENSIONAL_FACTOR.md) | f(1+f^2) = (d-1)/2 | 12 + 25 | VERIFIED |
| [D05](D05_DYNAMIC_MODE_COUNT.md) | Osc = 4^N-(N+1)-Stat(N) | 4 + 22 + 23 | VERIFIED |
| [D06](D06_SPECTRAL_GAP.md) | Gap = 2*gamma | 1 + 3 | VERIFIED |
| [D07](D07_Q_DISTRIBUTION.md) | Q ~ arcsine distribution | 2 + 7 | VERIFIED |
| [D08](D08_CROOKS_RATE_IDENTITY.md) | ln(d_fast/d_slow) = 2·artanh(Δd/(2Σγ)) | 1 | PROVEN |
| [D09](D09_SECTOR_SFF_PAIRING.md) | K_freq(w,t) = K_freq(N−w,t) | 1 | PROVEN |
| [D10](D10_W1_DISPERSION.md) | ω_k = 4J(1−cos(πk/N)) | Tight-binding reduction | PROVEN |

## Dependency Graph

```
Formula 2 (dispersion) ----> D1 (bandwidth)
                        \--> D7 (Q distribution)
                              |
Formula 7 (Q spectrum) ------/
                        \--> D2 (V = Q_max/Q_mean)
Formula 6 (V-Effect) -------/

Formula 27 (K per channel) -> D3 (crossing ratios)

Formula 12 (1-qubit cross) -> D4 (dimensional factor)
Formula 25 (Bell+ cross) --/

Formula 4 (Stat(N)) -------> D5 (oscillating count)
Formula 22 (XOR drain) ---/
Formula 23 (XOR vanish) -/

Formula 1 (palindrome) ----> D6 (spectral gap)
Formula 3 (decay bounds) -/
                        \--> D8 (Crooks rate identity)
                        \--> D9 (sector SFF pairing)

Tight-binding reduction -----> D10 (w=1 dispersion)
                                |---> Formula 7 (Q spectrum)
                                |---> Formula 41 (palindromic time)
                                \---> D1 (bandwidth)
```

## Failed Cascades

Explored but no closed formula found:

- **Cascade B** (Protection Factor): PF grows monotonically with
  contrast. No analytical optimum. Reason: Hamiltonian mixes
  weight-parity sectors (w with w+-2), making effective rates
  a matrix eigenvalue problem.

- **Cascade C** (Optimal Contrast): Depends on B. No universal
  closed form. Practical limit is hardware noise floor.

- **Cascade E** (2x Law vs Thermal): Amplitude damping changes
  the palindromic structure. Ratio depends on gamma_ad/gamma_deph.
  No universal threshold.

- **Cascade D** (Sweet Spot vs N): Trivial restatement of
  formulas 18 + D6. Q ~ 800*N at fold boundary.

## Verification Scripts

- [`simulations/verify_derivations.py`](../../../simulations/verify_derivations.py) (D1-D6)
- [`simulations/derivation_cascades.py`](../../../simulations/derivation_cascades.py) (D7 + failed cascades)
- Results: [`simulations/results/verify_derivations.txt`](../../../simulations/results/verify_derivations.txt)
- Results: [`simulations/results/derivation_cascades.txt`](../../../simulations/results/derivation_cascades.txt)
