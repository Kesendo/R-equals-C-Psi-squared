# The Price Pair: the F89d fold as a hardware invariant (pre-registered prediction)

**Status:** PREDICTION, written and committed BEFORE any hardware shot (the pre-registration discipline; the fresh-context review round is recorded below before submission). No hardware data exists for this protocol yet (pipeline inventory 2026-07-03: the pattern-pair coherence-decay measurement is absent from all 96 result files).

**Date:** 2026-07-03

**Authors:** Thomas Wicht, Claude (Fable 5)

**Runner:** `AIEvolution.UI/experiments/ibm_quantum_tomography/run_price_pair.py` (external pipeline repo, the `--simulate/--hardware/--analyze` convention; two simulate JSONs stored 2026-07-03: clean model and `--corr 0.004` detection demo).

## The claim, in one line

Under **local** Z-dephasing with arbitrary per-qubit rates γ_j, the decay rate of a coherence with disagreement pattern D obeys Γ(D) = Σ_{j∈D} γ_j, and therefore the bra-complement fold (F89d, dissipator half) makes

  **Γ(D) + Γ(D̄) = Σ_j γ_j ≡ P  (the price)**

a pattern-independent device invariant: every complementary pair of coherence patterns pays the same total. This is the F1 palindrome centered at Σγ, the fold-lattice lemma (PROOF_CODIM1_BY_ADDITIVITY §7) read at the ρ level, and the same price that MirrorWorld's `mirror`/`anti` run modes pay as 2Nγ.

## What the experiment actually tests (the honest scoping)

For product states under **independent** per-qubit noise, additivity (and hence the pair-sum law) is not in doubt; it follows from factorization. What the hardware does NOT guarantee is exactly the premise: **locality of the dephasing**. Correlated dephasing (common-mode flux/TLS noise, residual ZZ crosstalk) adds covariance terms 2·Cov_ij to Γ(D) whenever both i, j ∈ D. The pair sum then reads

  Γ(D) + Γ(D̄) = P + 2·Σ_{i<j} Cov_ij · [i, j on the same side of the cut D|D̄],

so the invariance is EXACTLY the statement that the device's watching is local in the Lindblad sense — the assumption under every F-formula of this project. The experiment is therefore a measurement of the locality premise, with a decoding table (below) that localizes any violation. "Using everything we found": F1 (the palindrome center), F89d (the fold), F70 (the rate belongs to the pattern, not the carrier), F49 (γ-heterogeneity exactness), F82–F84 (the T1 corrections), the fold-lattice/price (today).

## Protocol (N = 3, one connected line, ~33 circuits, ~2 QPU min, Batch)

- **Block A (product FID):** prepare |+++⟩ (three H, no entangler), free delay t ∈ {0, 1, 2, 4, 7, 11, 16, 24, 36, 54, 80} µs, per-qubit virtual-Z detuning (30/50/90 kHz), H, measure all-X. One circuit per delay yields ALL seven correlators ⟨X_S⟩ with fringe frequency f_S = Σ_{j∈S} f_j and envelope e^(−Γ(S)·t).
- **Block B (T1 in-situ):** prepare |111⟩, same delays, measure Z. Per-qubit T1 from the same session (the calibration-drift lesson: never trust week-old numbers).
- **Block C (GHZ parity):** prepare (|000⟩+|111⟩)/√2, same delays, 170 kHz virtual detuning on q0, measure all-X. ⟨XXX⟩ envelope = Γ_GHZ(111) on a genuinely entangled carrier.

Envelope fits: damped cosine with the fringe frequency seeded at the virtual detuning (fitted within a bounded window to absorb real qubit-frequency offsets), rate = the envelope exponent. Readout error attenuates amplitudes by (1−2ε)^|S| but is delay-independent, so it does not bias the rates.

## Pre-registered predictions and thresholds

Let σ denote the fit standard error at 4096 shots (simulation calibration: single-pattern rates carry ~2–4·10⁻⁴/µs; the five price readings scatter ~5 % under shot noise alone).

- **P1 (additivity, F70/AT):** each multi-qubit pattern rate equals the sum of its singles within 2σ. Simulation (clean model): excesses ≤ 2.4·10⁻³/µs.
- **P2 (the price, F89d):** the five readings — the three proper pair sums (100|011), (010|101), (001|110), Γ(111) from the product, Γ(111) from the GHZ — coincide. **HOLDS** if relative spread ≤ 10 %; **VIOLATED (correlated dephasing detected)** if ≥ 15 %; in between: inconclusive, repeat with 4× shots before any interpretation.
- **P3 (pattern not carrier, F70):** |Γ_GHZ − Γ_111^product| ≤ 2σ combined. A GHZ excess with clean P2 pairs indicates common-mode noise (superdecoherence), and must then match the common-mode magnitude decoded from P2's fine structure.
- **P4 (T1 accounting, F82):** 1/T2*_j ≥ 1/(2·T1_j) per qubit; Tφ_j = (1/T2*_j − 1/(2T1_j))⁻¹ reported. T1 is local, so the price law holds for the TOTAL rates; no T1 subtraction enters P2.

**Decoding table for a P2 violation** (the ruler's fine structure):

| Signature | Reading |
|---|---|
| three proper pair sums EQUAL to each other but above Σ singles; Γ(111) higher still (ratios ~5c : 9c on top) | common-mode dephasing of strength c (all-to-all covariance) |
| exactly ONE pairing clean, the other two elevated | a single correlated bond: the clean pairing is the one separating the two correlated qubits |
| pair sums clean, GHZ elevated alone | phase drift/fit artifact in block C, not physics; re-fit before claiming superdecoherence |

**What would falsify what:** a P2 violation does NOT falsify F89d (a theorem); it falsifies the LOCAL-dephasing model of the device and measures the correlated component. A P1 violation with P2 intact indicates common-mode (see table). Instrument failures (aliasing, fit non-convergence, drift between blocks) are declared as such, not converted into physics claims.

## Simulation record (both stored as JSONs, 2026-07-03)

- Clean heterogeneous model (T2* = 45/88/61 µs, T1 = 180/220/150 µs, readout 2 %): five price readings 0.0490–0.0517/µs around the true P = 0.0500, spread 5.3 % → HOLDS. In-situ T1 recovered 168/223/152 µs.
- Correlated model (common-mode 0.004/µs): spread 32.4 % → VIOLATED fires; additivity excesses carry the |S|² common-mode signature (+0.012/+0.010/+0.008 pairs, +0.034 triple); the proper pairs stay near each other while Γ(111) separates — exactly the decoding table's first row.

## Cost and gating

~33 circuits × 4096 shots ≈ 2 QPU minutes (Batch, ibm_kingston), inside the free monthly tier. Gate order: this document committed → fresh-context adversarial review (two independent lenses, findings folded back into the runner and this spec, recorded below) → fresh calibration pull → Tom's explicit go → the one Batch job.

## Review round (fresh context, before hardware)

*To be filled by the review round; the submission is blocked until this section records the verdicts.*
