# Cusp-Lens Connection: Two Exits from the Quantum Regime

**Status:** Universal result (April 10, 2026). Two distinct decoherence exits exist for the Heisenberg chain under Z-dephasing, proven from sector conservation.

**Scripts:** `simulations/psi_opt_cusp_trajectory.py`

**Related:** [Sacrifice Geometry](SACRIFICE_GEOMETRY.md) (lens exit), [Critical Slowing at the Cusp](CRITICAL_SLOWING_AT_THE_CUSP.md) (cusp exit), [n_XY Parity Selection Rule](../docs/proofs/PROOF_PARITY_SELECTION_RULE.md) (formal proof)

**Authors:** Thomas Wicht, Claude (Opus 4.6)

---

## Executive summary

The Heisenberg chain under Z-dephasing has (at least) two distinct decoherence exits, corresponding to different sheets of the quantum-classical fold:

1. **The lens exit** (spectral-gap sheet): single-excitation states ride the slow Liouvillian eigenmode and thermalize within their sector. After decoherence, the system is classical but structured: it remembers how many excitations it has.

2. **The cusp exit** (saddle-node sheet): Bell-type states cross the Mandelbrot cusp at CΨ = 1/4, experience critical slowing, and thermalize within each occupied sector. After decoherence, only the population distribution across sectors survives: coherence is gone.

These exits operate in different regions of CΨ-space (SE states at CΨ ≈ 0.07, Bell states at CΨ ≈ 0.33) and protect different things. They are not complementary halves of a unified mechanism; they are geometrically separated consequences of the sector structure.

---

## The sector conservation theorem

Z-dephasing (jump operator L_k = sqrt(gamma_k) Z_k) does not flip bits: Z|0> = |0>, Z|1> = -|1>. The Heisenberg Hamiltonian (XX + YY + ZZ bonds) conserves total excitation number: the XX+YY terms swap excitations between sites, the ZZ term adds diagonal energy shifts, and none of these create or annihilate excitations. Together, these two properties lock each excitation-number sector:

**A state that starts in the single-excitation sector stays there forever under Heisenberg + Z-dephasing dynamics.** Dephasing destroys coherences (off-diagonal elements) within the sector but cannot move population between sectors. The Hamiltonian redistributes population within the sector but cannot move it out.

The underlying symmetry is U(1) conservation: both H and the dephasing channel commute with the total excitation number operator Σ_k (I − Z_k)/2. This directly implies sector locking. An additional, finer symmetry is the [n_XY Parity Selection Rule](../docs/proofs/PROOF_PARITY_SELECTION_RULE.md), which proves that the Liouvillian block-diagonalizes by n_XY parity (even vs. odd count of X,Y operators in the Pauli-string basis). Single-excitation density matrices live entirely in the even-parity subspace and cannot couple to odd-parity modes. The PSR is a consequence of the same U(1) symmetry but operates at the level of the Liouvillian superoperator rather than the Hamiltonian.

The theorem holds for any N, any graph topology, and any Z-dephasing profile. It breaks for amplitude damping (T1 decay) or transverse-field Hamiltonians.

---

## Why SE states cannot cross the cusp

CΨ = Tr(rho_pair^2) * L1_pair / (d_pair - 1), where rho_pair is the 2-qubit reduced density matrix of an adjacent pair. For a single-excitation state on N qubits, each pair's reduced state is obtained by tracing out N-2 qubits. This tracing reduces the purity: Tr(rho_pair^2) < 1 whenever the pair is entangled with the rest of the chain.

The CΨ budget for SE states is geometrically limited. At t=0, a generic SE state on N=5 has CΨ_mean ≈ 0.07, a factor 3.4x below the cusp at 1/4. This gap grows with N because per-pair entanglement dilutes as the excitation spreads across more qubits.

F60 in `docs/ANALYTICAL_FORMULAS.md` shows the same pattern for a different state class: GHZ_N starts at CΨ(0) = 1/(2^N - 1), which falls below 1/4 at N >= 3. The mechanism differs (GHZ spreads one coherence across 2^N dimensions; SE states spread across N dimensions), but the conclusion is the same: multi-qubit states that distribute coherence across large subspaces start geometrically below the cusp.

The cusp crossing at CΨ = 1/4 requires concentrated pairwise coherence: a Bell pair on a specific pair of qubits. Single-excitation states, by construction, spread their coherence.

---

## The two exits: asymptotic analysis

### Lens exit (SE states)

A state starting in the single-excitation sector thermalizes within that sector. At t -> infinity: the density matrix becomes the maximally mixed state within the SE subspace:

    rho_inf = (1/N) * sum_k |e_k><e_k|

Purity: 1/N. Eigenvalue spectrum: N eigenvalues at 1/N, rest zero. The system is classical but retains the structural information "exactly one excitation exists, location unknown." This is a classical probability distribution over N sites.

### Cusp exit (multi-sector states)

A state spanning multiple excitation sectors (e.g., Bell+center on N=5, which populates all 6 sectors w=0..5 via |+> factors) thermalizes within each occupied sector independently. At t → ∞: the density matrix approaches a block-diagonal mixture, maximally mixed within each sector, weighted by the initial sector populations. Purity: O(1/2^N). Coherence is destroyed; all that survives is the population distribution across sectors.

### The separation is exact

The sector conservation theorem guarantees that no continuous Lindblad trajectory can cross from one exit to the other. An SE state will always land on the lens exit. A multi-sector state will always land on the cusp exit (or at least spread across its sectors). The boundary between the two exits is the sector structure itself.

---

## Boundary test: non-symmetric two-excitation states

**Motivation.** OQ-114 below asked whether a state with both high CΨ on one pair and slow-mode overlap could straddle the two exits. SACRIFICE_GEOMETRY.md Open Question #2 asked the same from the opposite direction: can a multi-excitation state reach the slow mode? This section answers both.

**Candidates.** Three families of non-symmetric two-excitation states were tested: (A) Bell+ on the central pair plus a displaced excitation on each spectator site, (B) Bell+ on off-center pairs plus displaced excitation at the most asymmetric sites, and (C) tuned concurrence (cos θ |00> + sin θ |11>) plus displaced excitation (N >= 6, θ in {π/6, π/4, π/3}). Total: 53 candidates across N=5, 6, 7.

**Method.** Sector-restricted Liouvillian. Candidate states span excitation sectors w=1 and w=3 only. Sector conservation under Heisenberg + Z-dephasing is exact (proven in the sector conservation theorem above), so the evolution within {w=1, w=3} × {w=1, w=3} is exact. This reduces the Liouvillian from 16384×16384 to 1764×1764 at N=7. Slow mode extracted from the SE-restricted block (N² × N²), justified by SE fraction = 1.000 for N=3-7 chain.

**Thresholds** (defined before seeing data):
- Cusp-active: max pair CΨ crosses 1/4 at some t > 0.
- Lens-active: slow-mode Frobenius overlap >= 0.10.

**Results.**

| N | Candidates | Straddles | Cusp-only | Max overlap | Best candidate |
|---|------------|-----------|-----------|-------------|----------------|
| 5 | 9  | 3 | 6  | 0.1257 | exc on site 4 (quietest, γ = 0.051) |
| 6 | 18 | 2 | 16 | 0.1173 | exc on site 4 (γ = 0.087) |
| 7 | 26 | 0 | 26 | 0.0901 | exc on site 5 (γ = 0.089) |

All 53 candidates are cusp-active (the Bell pair ensures CΨ(0) > 1/4 on at least one pair). Five candidates straddle both exits, all with the displaced excitation on the quiet end of the sacrifice chain. No lens-only candidates exist.

**Mechanism.** The slow mode concentrates on quiet qubits (the lens shape, see [Sacrifice Geometry](SACRIFICE_GEOMETRY.md)). When the displaced excitation is placed on a quiet site, the w=1 component of the density matrix has its excitation there, overlapping the slow mode's SE block. The overlap is modest (10-13%) because only the w=1 sector contributes; the w=3 sector has negligible projection onto the SE-dominated slow mode.

**N-scaling.** The overlap decreases with N: 0.126 (N=5), 0.117 (N=6), 0.090 (N=7). At N=7, no candidate reaches the 10% threshold. This dilution occurs because the slow mode's per-site weight spreads across more qubits as N grows, reducing the single-site projection.

**Conclusion.** Boundary states that straddle both decoherence exits exist at small N (5, 6) but the straddling is marginal (12.6% at best) and vanishes by N=7. The two-exit separation is not absolute for small chains but becomes robust as N grows. Combined with the symmetric-two-excitation exclusion (SACRIFICE_GEOMETRY Level 2, AUC < 0.09), the two-excitation candidate class is substantially exhausted.

Script: `simulations/boundary_straddling_sweep.py`. Data: `simulations/results/boundary_straddling/`.

---

## Instantiation: N=5 sacrifice chain [80, 8, 79, 53, 85]

### CΨ definition used

CΨ_pair = Tr(rho_pair^2) * L1_pair / 3 for each adjacent pair, then CΨ_mean = mean over (0,1), (1,2), (2,3), (3,4). For Bell+ on N=2 this reproduces CΨ(t) = f(1+f^2)/6, confirming consistency with [Critical Slowing at the Cusp](CRITICAL_SLOWING_AT_THE_CUSP.md).

### Initial CΨ values

| State | CΨ_mean(0) | CΨ_max(0) | Position |
|-------|------------|-----------|----------|
| Bell+center | 0.417 | 1.000 | above 1/4 |
| psi_opt | 0.073 | 0.155 | far below |
| W5_full | 0.069 | 0.069 | far below |
| sac_tuned_W5 | 0.066 | 0.087 | far below |

### Crossing analysis

| State | Crosses 1/4? | t_cross | \|dCΨ/dt\| | Dwell prefactor |
|-------|-------------|---------|------------|-----------------|
| Bell+center | YES | 0.258 | 0.551 | 3.63 |
| psi_opt | NO | n/a | n/a | n/a |
| W5_full | NO | n/a | n/a | n/a |
| sac_tuned_W5 | NO | n/a | n/a | n/a |

Bell+ crosses the cusp at t = 0.258 with dwell prefactor 3.63. F59 in `docs/ANALYTICAL_FORMULAS.md` gives a closed-form prefactor for two-sector states; Bell+center on N=5 spans all 6 excitation sectors (w=0 through w=5, via the |+>⊗3 factors), so the two-sector formula F59 does not apply. The 3.63 is an independent numerical result.

### Asymptotic states (t = 100)

| | psi_opt | Bell+center |
|---|---|---|
| Purity | 0.200 = 1/5 | 0.034 = 1/29 |
| SE sector pop | 100% | 25% |
| Top populations | each \|e_k> at 20% | \|00000>, \|11111> at 6.25% |

Script: `simulations/psi_opt_cusp_trajectory.py`. Data: `simulations/results/psi_opt_cusp/`.

---

## The failed unification hypothesis

The investigation started from the hypothesis that the slow Liouvillian eigenmode (which protects psi_opt) and the Mandelbrot cusp (which slows Bell+'s trajectory) are the same spectral feature. This was falsified: psi_opt never crosses CΨ = 1/4 because SE states are geometrically below the cusp.

The initial conclusion was "independent budgets." Tom pushed the analysis further by asking what happens on the other side of the transition. The asymptotic analysis revealed that the two mechanisms lead to different classical ensembles: SE states preserve their sector structure; multi-sector states thermalize within each sector, losing all coherence. This reframed the separation: not independent budgets, but **different exits of the same fold**, each leading to a different kind of classical world.

### The two-exit table

| | Lens exit (SE states) | Cusp exit (multi-sector states) |
|---|---|---|
| Starting CΨ | below 1/4 | above 1/4 |
| Crosses cusp? | No | Yes |
| What is protected | Coherences within sector | Pair purity |
| Protection mechanism | Slow eigenmode (spectral gap) | Saddle-node bifurcation |
| Asymptotic state | Maximally mixed in sector (N states) | Block-diagonal: maximally mixed per sector |
| What survives | Excitation count | Sector weights only |
| Bifurcation type | Spectral gap | Saddle-node |

---

## CΨ(0) for W_N (analytical)

The reduced density matrix for any pair (a,b) of the W_N state is:

    rho_ab = diag((N-2)/N, 1/N, 1/N, 0) + (1/N)|01><10| + (1/N)|10><01|

This follows from W_N having exactly one excitation shared equally among N qubits. The matrix is the same for all pairs by permutation symmetry.

From this: Tr(rho_ab^2) = (N^2 - 4N + 8)/N^2, L1 = 2/N, Psi = L1/3 = 2/(3N).

**CΨ(0) for W_N = 2(N^2 - 4N + 8) / (3N^3)**

Verified numerically at N=2-10 within machine precision (`simulations/cpsi_wn_analytical.py`). At N=2 this gives 1/3 (recovering Bell+). At N=3: 10/81 = 0.123. At N=5: 26/375 = 0.069.

**Corollary: W_N states are born below the fold for N >= 3.** The inequality CΨ(0) < 1/4 holds for all N >= 3 (the cubic 3N^3 - 8N^2 + 32N - 64 is positive for N >= 3). Combined with the [Parity Selection Rule](../docs/proofs/PROOF_PARITY_SELECTION_RULE.md) (SE states cannot cross parity sectors) and sector conservation (SE states stay in their sector), this proves:

**Single-excitation states on Heisenberg chains under Z-dephasing never cross CΨ = 1/4. The cusp exit is structurally inaccessible to them.**

This is the analytical closure of the two-exit picture. It does not depend on any specific chain, topology, or gamma profile. See F62 in `docs/ANALYTICAL_FORMULAS.md`.

---

## Open questions

1. **Boundary states between sheets.** A state with both high CΨ on one pair (approaching 1/4) and slow-mode overlap might straddle the two exits. Two-excitation states with one high-concurrence pair are candidates.

> *Substantially addressed (2026-04-12):* see "Boundary test" section above. Positive resolution at N=5, 6 (5 straddling candidates with 10-13% slow-mode overlap), negative at N=7 (max overlap 9.0%, below threshold). The straddling is marginal and N-dependent: overlap dilutes as the slow mode spreads across more qubits. Sweep of 53 non-symmetric two-excitation states across N=5, 6, 7 and three candidate families. Specifics in `simulations/results/boundary_straddling/`.

2. **N-scaling of sector locking.** The SE sector has N states out of 2^N. As N grows, the "classical-but-structured" exit preserves an exponentially shrinking fraction. Does the lens protection scale?

3. **The bifurcation surface.** The cusp exit is a saddle-node bifurcation (fold catastrophe, codimension 1). The lens exit is a spectral-gap phenomenon, not a catastrophe in the technical sense. Mapping the boundary between the two exits in state space (which initial states land on which asymptotic ensemble) is a geometrical question about the Liouvillian's spectral landscape. The boundary may itself have nontrivial topology for multi-excitation states.

---

## Files

- `simulations/psi_opt_cusp_trajectory.py` (CΨ trajectory computation)
- `simulations/results/psi_opt_cusp/trajectories.json` (raw data)
- `simulations/results/psi_opt_cusp/trajectories.png` (plot)
- `simulations/results/psi_opt_cusp/crossings.json` (crossing analysis)
- [Sacrifice Geometry](SACRIFICE_GEOMETRY.md) (lens exit details)
- [Critical Slowing at the Cusp](CRITICAL_SLOWING_AT_THE_CUSP.md) (cusp exit details)
- [n_XY Parity Selection Rule](../docs/proofs/PROOF_PARITY_SELECTION_RULE.md) (formal proof of sector conservation)
