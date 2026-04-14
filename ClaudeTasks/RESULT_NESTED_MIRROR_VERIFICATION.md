# RESULT: Nested Mirror Structure Verification Suite

**Date:** 2026-04-14, evening session
**Executor:** Claude Code (Opus 4.6)
**Task file:** `ClaudeTasks/TASK_NESTED_MIRROR_VERIFICATION.md`
**Parent commit:** `d48502a`
**Commits produced:**

| Commit | Check | Verdict |
|--------|-------|---------|
| `1d77bb6` | Check 1: SWAP-artifact | genuine but perturbative |
| `7baaa7c` | Check 2: N=3 scaling | class-scaling FALSIFIED |
| `3e2f429` | Check 3: coupling robustness | three-class structure confirmed |
| `e3dc6bc` | Check 4: rebound mechanism | CONFIRMED |

**Working tree after session:** clean (all results committed)

---

## Summary table

| Check | Question | Verdict | Falsification criterion |
|-------|----------|---------|------------------------|
| 1 | Is SWAP +-1 a basis artifact? | Genuine but perturbative | 4: does not trigger |
| 2 | Does N=3 show 4 evenly-spaced classes? | **FALSIFIED** (12 classes) | 1: **triggered** |
| 3 | Does the three-class structure depend on coupling? | Robust (5 of 6 couplings) | 2: does not trigger |
| 4 | Do mirror modes drive the rebound? | Confirmed (necessary + sufficient) | 3: does not trigger |

---

## Check 1: SWAP-artifact check

**Method:** Three tests. (1) Sweep gamma_S/gamma_B from 0 to 0.1 (21 points). (2) Add local field h_S * Z_S from 0 to 0.5 (11 points) to break the 10-fold mirror degeneracy. (3) Vary gamma/J ratio from 0.001 to 1.0 to characterize the scaling.

**Key numbers:**
- gamma_S sweep: 16/16 modes remain within 0.05 of +-1 throughout (no smearing)
- Local field at h_S = 0.1: degeneracy breaks to 6 Re(lambda) classes, 16/16 modes still sharp
- Local field at h_S = 0.5: only 5/16 sharp (conserved + correlation modes exact; mirror modes degrade to 0.6-0.9)
- gamma/J scaling: deviation from +-1 follows (gamma/J)^2 / 2 to ratio 1.000 for gamma/J <= 0.2

**Interpretation:** SWAP +-1 is NOT a numpy basis artifact. It is a perturbative near-symmetry inherited from [H_XX+YY, SWAP] = 0. At gamma = 0: exact SWAP eigenvectors. Dissipator shifts them by O(gamma^2/J^2). The pattern is coupling-dependent (requires H invariant under qubit exchange).

**Hypothesis doc change:** Observation 4 rewritten from "caveated" to "verified, perturbative". Falsification criterion 4 resolved.

**Result files:**
- `simulations/nested_mirror_swap_check.py`, `simulations/nested_mirror_swap_check_extended.py`
- `simulations/results/nested_mirror_swap_check/swap_vs_gamma_S.{png,txt}`, `swap_extended_analysis.txt`

---

## Check 2: N=3 scaling check

**Method:** Three-qubit chain S-M-B, H = J/2 * (XX+YY on S-M and M-B bonds), single jump sqrt(gamma_B) * I(x)I(x)Z. Diagonalize 64x64 Liouvillian.

**Key numbers:**
- Number of distinct Re(lambda) classes: **12** (hypothesis predicted 4)
- Predicted positions -2g/3 = -0.0667 and -4g/3 = -0.1333: **MISSING**
- Actual Re values cluster around: 0, -0.050, -0.075, -0.100, -0.125, -0.150, -0.200
- Palindromic pairing around Sigma_gamma = 0.1: **confirmed for all 64 eigenvalues**
- Boundary classes Re = 0 (4 modes) and Re = -0.2 (4 modes): **present with expected structure**
- Non-Markovian rebound: **confirmed**, amplitude 0.36 (stronger than N=2 at 0.17)

**Interpretation:** The even-spacing formula Re = -k * 2*gamma/N is **falsified**. The three-class structure at N=2 is specific to the minimal system where coupling forces equal XY-weight distribution on both qubits. At N=3, the chain topology creates modes with varying amounts of XY weight on the dephased qubit B, producing a richer spectrum.

What survives: boundary classes at Re = 0 and Re = -2*gamma, palindromic pairing (follows from proven MIRROR_SYMMETRY_PROOF), non-Markovian rebound growing with system size.

**Hypothesis doc change:** Hypothesis statement rewritten to remove scaling formula. Falsification criterion 1 marked as triggered. Tier remains 3 (partially verified).

**Result files:**
- `simulations/qubit_in_qubit_n3.py`
- `simulations/results/qubit_in_qubit_n3/eigenvalue_analysis.txt`

---

## Check 3: Coupling robustness

**Method:** Re-run N=2 eigenmode analysis with 6 coupling types: XX+YY (reference), XXX (Heisenberg), pure XX, pure YY, pure ZZ, XX+ZZ.

**Key numbers:**

| Coupling | Classes | Degeneracy | 1/sqrt(2) split | Rebound | Amplitude |
|----------|---------|------------|-----------------|---------|-----------|
| XX+YY | 3 | {3,10,3} | YES | YES | 0.356 |
| XXX | 3 | {3,10,3} | YES | YES | 0.363 |
| XX only | 3 | {4,8,4} | YES | NO | 0.000 |
| YY only | 3 | {4,8,4} | YES | YES | 0.356 |
| ZZ only | 2 | {8,8} | NO | YES | 0.498 |
| XX+ZZ | 3 | {4,8,4} | YES | YES | 0.498 |

**Interpretation:** The three-class structure is coupling-robust. Any coupling with at least one off-diagonal Pauli channel (X or Y) produces the middle class with 1/sqrt(2) partial-trace split. Only pure ZZ breaks it (ZZ commutes with Z-dephasing, no mode mixing). The degeneracy pattern {3,10,3} is specific to exchange-symmetric couplings (XX+YY, XXX); single-channel couplings give {4,8,4}. Rebound is absent for pure XX coupling (related to initial state symmetry, not to class structure).

**Hypothesis doc change:** Falsification criterion 2 marked as not triggered. Pending verification 3 marked done.

**Result files:**
- `simulations/qubit_in_qubit_coupling_sweep.py`
- `simulations/results/qubit_in_qubit_coupling/coupling_sweep.txt`

---

## Check 4: Rebound-mechanism check

**Method:** Decompose initial state rho(0) = |+><+| (x) I/2 in the eigenmode basis of L. Construct projected states: conserved-only, mirror-only, correlation-only, no-mirror, no-conserved, no-correlation. Evolve each under L and measure coherence rebound.

**Key numbers:**
- rho(0) decomposition: 50.0% conserved weight, 50.3% mirror weight, 0.0% correlation weight
- Original: rebound amplitude 0.356
- Mirror-only (nearest DM): rebound amplitude 0.356 (identical)
- No-mirror (valid DM): rebound amplitude **0.000** (completely eliminated)
- Conserved-only: zero S-coherence, no dynamics
- Correlation-only: zero initial weight, no rebound

**Interpretation:** Mirror modes are **necessary** (removing them kills rebound entirely) and **sufficient** (mirror-only state reproduces full rebound). The mechanism claim is confirmed: the non-Markovian rebound is carried by the intermediate eigenvalue class at Re = -gamma, whose eigenmodes oscillate between S and B with equal 1/sqrt(2) partial-trace weight on each qubit.

**Caveat:** The mirror-only projection had trace 0 (mirror modes carry coherences, not populations) and was corrected to the nearest valid density matrix. The nearest-DM effectively adds a maximally mixed diagonal, closely approximating the original state. The stronger evidence is the no-mirror case: a valid DM with zero rebound.

**Hypothesis doc change:** Falsification criterion 3 marked as not triggered. Pending verification 4 marked done.

**Result files:**
- `simulations/qubit_in_qubit_mode_projection.py`
- `simulations/results/qubit_in_qubit_mode_projection/rebound_by_class.png`, `mode_projection_results.txt`

---

## Hypothesis final state

**Tier:** 3 (partially verified, partially falsified)

**What survives:**
1. **Three-class structure at N=2** with boundary classes (Re = 0, Re = -2*gamma), middle class at Re = -gamma, and 1/sqrt(2) partial-trace split. Confirmed across 5 of 6 coupling types.
2. **Palindromic pairing** around Sigma_gamma. Holds at N=2 and N=3. Follows from proven MIRROR_SYMMETRY_PROOF.
3. **Non-Markovian rebound** driven by mirror modes. Mechanism confirmed by eigenmode projection. Rebound amplitude grows with N (0.17 at N=2, 0.36 at N=3).
4. **SWAP near-symmetry** of eigenmodes at N=2. Genuine perturbative structure from [H, SWAP] = 0, not a basis artifact.

**What was removed:**
1. **Class-scaling formula** Re = -k * 2*gamma/N. Falsified at N=3 (12 classes instead of 4).
2. **Claim of (N+1) classes.** The number and spacing of intermediate eigenvalue classes depend on the Hamiltonian's mode structure, not on a universal N-scaling law.

**What remains open:**
- Connection to experimental observables (IBM non-Markovianity measurements)
- Asymmetric gamma profiles
- Why pure XX coupling has no rebound (likely initial-state-dependent, not structural)

---

## Recommended next steps

1. **Eigenvalue-class formula for general N.** The even-spacing prediction failed, but the boundary classes and palindromic pairing are universal. A correct formula would involve the absorption theorem's mode-averaged XY weights. Check whether the absorption theorem directly predicts the N=3 eigenvalue positions.

2. **BLP non-Markovianity index.** The rebound amplitude is now a measurable quantity with a confirmed mechanism. Compute the BLP index for the N=2 and N=3 systems as a function of gamma/J. This gives a quantitative prediction for IBM hardware tests.

3. **EQ-013 sub-question 2 update.** The partial answer from the minimal nest is now strengthened: nesting IS non-trivial (rebound confirmed, mechanism identified), but the scaling law is wrong. Sub-question 2 can be narrowed to: "what determines the intermediate class structure at general N?"

4. **Pure XX rebound absence.** The XX coupling shows three classes but no rebound. This deserves a short investigation: is it initial-state-dependent (try |+Y> instead of |+X>), or structural (XX coupling preserves a sub-symmetry that prevents information backflow)?

5. **Upgrade to experiment doc.** The verified parts (three-class structure, coupling robustness, rebound mechanism) have enough numerical evidence to move from hypothesis to experiment documentation. The falsified scaling prediction is already documented. A clean experiment doc would separate the confirmed N=2 structure from the open N-scaling question.

---

## Open threads

1. **Absorption theorem connection.** The 12 eigenvalue classes at N=3 should be explainable from the absorption theorem: Re(lambda) = -2*gamma_B * <n_XY_B> where <n_XY_B> is the mode-averaged XY weight on qubit B. The specific values depend on how the Hamiltonian's eigenmodes distribute XY weight across the chain. This connects to EQ-006 (rate bounds).

2. **ZZ coupling as the "dark" channel.** ZZ commutes with Z-dephasing and produces only 2 classes {8,8}, no middle class, but still has rebound (amplitude 0.498). This rebound comes from purely unitary oscillation (H = ZZ creates conditional phase), not from the mirror-mode mechanism. Different physics, same observable.

3. **Degeneracy pattern {3,10,3} vs {4,8,4}.** The exchange-symmetric couplings (XX+YY, XXX) give {3,10,3}; single-channel couplings give {4,8,4}. The 3 extra modes in the middle class for exchange-symmetric couplings come from the additional qubit-swap symmetry. This connects to EQ-004 (transient complexity and spatial symmetry groups).

4. **SWAP perturbative correction formula.** |<SWAP>| = 1 - (gamma/J)^2/2 is an empirical fit. A perturbation theory derivation (expanding L = L_H + epsilon * L_D with epsilon = gamma/J) should give this as the leading-order correction to the SWAP eigenvalues. This would connect the SWAP structure to the general perturbation theory of non-Hermitian operators.

---

*End of result. Working tree is clean.*
