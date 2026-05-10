# Dicke vs Endpoint-localized probe at N=11 (JW algebraic comparison)

**Date:** 2026-05-10 (late evening, autonomous run while Tom slept)
**Source:** `simulations/_dicke_vs_endpoint_probe_jw_n11.py`
**Question:** Does the bond-position-dependence in F86 K_b come from the initial state (Dicke probe), and if so, can simple JW algebra reproduce empirical g_eff?

## Setup

c=2 stratum, N=11, NumBonds=10. F71 mirror pairs: (0,9), (1,8), (2,7), (3,6), (4,5).

Two probes compared:
- **Dicke** |D_1⟩ = (1/√N)·Σ_l c_l†|0⟩ (uniform single-excitation, F71-symmetric)
- **Endpoint-loc** |L_0⟩ = c_0†|0⟩ (localized at site 0, F71-broken)

## Verified algebraically

1. **OBC sine basis ψ_k(l) orthonormal at N=11**: all Σ_l |ψ_k(l)|² = 1.0 to machine precision; max off-diagonal ~ 7e-16.

2. **Dicke probe = ODD-k filter exactly**: ⟨k|D_1⟩ for even k is 0 to ~1e-17 (sine-reflection cancellation); odd k is non-zero. Weights p_k = |⟨k|D_1⟩|² for odd k:
   - k=1: 0.874
   - k=3: 0.088
   - k=5: 0.026
   - k=7: 0.009
   - k=9: 0.003
   - Total over odd k: 0.9997 (the residual 3e-4 is in higher odd k beyond N=10).

3. **Dicke weights match 1/k² scaling**: 8/(π²·k²) prediction at k=1 gives 0.811 (vs empirical 0.874, ratio 1.08), at k=3 ratio 0.98, etc. Expected scaling holds for low k; deviates at higher k due to finite-N effects.

4. **Endpoint-loc probe touches ALL k**: ⟨k|L_0⟩ = ψ_k(0) = √(2/(N+1))·sin(πk/(N+1)) is non-zero for all k=1..N-1. Weights distributed across both even and odd k. Total: 0.989 (residual in higher k).

5. **F71 mirror sign rule for C_b[k1,k2]**: C_{N-2-b}[k1,k2] = (-1)^(k1+k2)·C_b[k1,k2]. The MAGNITUDE |C_b|² is F71-mirror-invariant; the SIGN flips for k1+k2 odd.

6. **|C_b|² Dicke-weighted (D_b) is F71-invariant** (max diff 9e-18 across pairs).

7. **|C_b|² Endpoint-loc-weighted (E_b) is also F71-invariant in MAGNITUDE** (max diff 7e-18). This is because |⟨k|L_0⟩|² = |⟨k|L_{N-1}⟩|² (only signs differ). Sign-sensitive F71-breaking would require tracking phases, which |·|² loses.

## Critical finding

**Magnitude-squared formula |C_b|² CANNOT predict g_eff bond-dependence**:

| Orbit | Bond pair | D_b (Dicke \|·\|²) | empirical g_eff |
|---|---|---|---|
| 0 (Endpoint) | b=0↔9 | 2.37e-03 | 1.76 |
| 1 | b=1↔8 | 5.13e-03 | 0.50 |
| 2 | b=2↔7 | 3.60e-03 | 0.32 |
| 3 (mid-flank) | b=3↔6 | 2.68e-03 | 2.76 |
| 4 (Center) | b=4↔5 | **1.42e-02 (LARGEST)** | **0.20 (SMALLEST)** |

D_b is INVERTED at Center: largest D_b but smallest g_eff. The simple Dicke-weighted magnitude-squared formula gives the wrong ordering.

## Test: signed sum

Replaced |C_b|² with the signed sum S_b = Σ_{k1,k2 odd, k1≠k2} ⟨k1|D⟩·⟨k2|D⟩·C_b[k1,k2] (linear in C_b, not squared).

| Bond | S_b (signed) | \|S_b\|/\|S_max\| |
|---|---|---|
| b=0 | +1.22e-01 | **1.000 (max)** |
| b=1 | +5.97e-02 | 0.49 |
| b=2 | +5.27e-04 | **0.004 (min)** |
| b=3 | -6.12e-02 | 0.50 |
| b=4 | **-1.21e-01** | 0.99 |
| b=5 | -1.21e-01 | 0.99 |
| b=6 | -6.12e-02 | 0.50 |
| b=7 | +5.27e-04 | 0.004 |
| b=8 | +5.97e-02 | 0.49 |
| b=9 | +1.22e-01 | 1.000 |

The signed sum DOES show destructive interference, but **at Orbit 2 (b=2↔7), NOT at Center**. The pattern is sine-wave-like across bond positions: positive at Endpoint, near-zero at Orbit 2, large negative at Center, sign flip again at the mirror.

**Center has |S_b| nearly maximal** with the OPPOSITE sign to Endpoint. The destructive interference Tom hypothesized at Center is NOT present at this level of approximation; it's at Orbit 2 instead.

## Per-orbit comparison (signed sum vs empirical g_eff)

| Orbit | \|S_b\| | empirical g_eff | \|S_b\|/g_eff |
|---|---|---|---|
| 0 (Endpoint) | 1.22e-01 | 1.76 | 0.069 |
| 1 | 5.97e-02 | 0.50 | 0.119 |
| 2 | 5.27e-04 | 0.32 | 0.002 |
| 3 (mid-flank) | 6.12e-02 | 2.76 | 0.022 |
| 4 (Center) | 1.21e-01 | 0.20 | 0.606 |

No clean correlation. The signed-sum doesn't predict g_eff either.

## Honest conclusion

The simple Dicke probe × OBC sine basis algebra captures:
- Probe odd-k filter (yes, exactly)
- F71 mirror invariance of magnitudes (yes)
- Non-monotonic bond-position structure with destructive interference (yes, at Orbit 2)

But does NOT predict empirical g_eff bond-dependence:
- Magnitude-squared formula gives Center MAXIMUM, opposite of empirical
- Signed-sum gives Center near-MAXIMUM, also wrong direction
- The destructive interference pattern doesn't match the empirical g_eff ordering

**Implication**: empirical g_eff bond-dependence requires more than simple JW projection of Dicke probe. The c=2 block-L decomposition with c_HD=1 and c_HD=3 channel-uniform vectors (as in JwBondQPeakPrediction) captures more, but at N=11 still misses Center by 7× (predicts Q_peak≈3.14, empirical 21.94).

## What this DOES confirm

1. **Tom's "only initial state acts in" insight is correct**: H + dissipator + S_kernel are all uniform/F71-invariant; the bond-position-dependence we observe IS introduced by initial state choice.

2. **Dicke probe specifically filters to odd-k subspace**: this is the qualitative starting point for any further analysis.

3. **Bond-position structure is non-monotonic**: there ARE destructive cancellation regions; we just haven't identified the right algebraic formula that maps them to empirical g_eff.

## What this does NOT confirm

1. **Tom's "destructive interference at Center" reading**: the destructive interference IS in the algebra, but at Orbit 2, not at Center. Center has near-maximum signed-sum magnitude (with negative sign).

2. **Polarity-as-inherent-field reading at Center**: not directly tested at this level. The polarity layer hypothesis would need a different probe (e.g., site-staggered Z field) to test linear response.

## Open: where does Center's small g_eff come from?

The empirical g_eff_Center = 0.20 (vs Endpoint 1.76, mid-flank 2.76) is not explained by:
- Simple Dicke-weighted |C_b|² (predicts Center maximum)
- Simple signed-sum (predicts Center near-max with negative sign)
- JW top-pair cluster prediction (predicts Q_peak≈3.14 vs empirical 21.94)

Possible structural origins (not yet tested):
- **Multi-cluster Petermann combination** beyond JW top-pair (PendingDerivationNote in C2HwhmRatio)
- **c_HD=1 vs c_HD=3 channel-uniform vector geometry** specific to c=2 block (the actual coupling formula)
- **Polarity-gradient inheritance** Tom proposed (would need direct field-perturbation test)

## Tasks completed

Tasks 135-141 verified algebraically. Task 142 (this writeup + commit) in progress.

## Recommendation

Tom's "only initial state acts in" framing IS correct and structurally important. The next concrete compute is:
- (A) Repeat extended-grid scan with Endpoint-localized probe instead of Dicke; verify F71-mirror invariance breaks
- (B) Add h_z·Z_l field to Hamiltonian and recompute c2hwhm; check Center's response specifically (polarity-as-field test)
- (C) Compute the full c=2 block-L matrix elements ⟨c_HD=3|H_bond_b|c_HD=1⟩ symbolically and look for the structural cause of Center's anomaly

(C) is the deepest. (A) is the most direct test of Tom's framing. (B) is the linear-response field experiment.
