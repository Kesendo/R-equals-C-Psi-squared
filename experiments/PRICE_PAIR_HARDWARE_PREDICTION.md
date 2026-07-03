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

so the invariance is, within the Gaussian-dephasing covariance model, equivalent to the device's watching being local in the Lindblad sense (and outside that model sensitive to further non-localities with their own signatures, see the decoding table): the assumption under every F-formula of this project. The experiment is therefore a measurement of the locality premise, with a decoding table (below) that localizes any violation. "Using everything we found": F1 (the palindrome center), F89d (the fold), F70 (the rate belongs to the pattern, not the carrier), F49 (γ-heterogeneity exactness), F82–F84 (the T1 corrections), the fold-lattice/price (today).

## Protocol (N = 3, one connected line, ~33 circuits, ~2 QPU min, Batch)

- **Block A (product FID):** prepare |+++⟩ (three H, no entangler), free delay t ∈ {0, 1, 2, 4, 7, 11, 16, 24, 36, 54, 80} µs, per-qubit virtual-Z detuning (30/50/90 kHz), H, measure all-X. One circuit per delay yields ALL seven correlators ⟨X_S⟩ with fringe frequency f_S = Σ_{j∈S} f_j and envelope e^(−Γ(S)·t).
- **Block B (T1 in-situ):** prepare |111⟩, same delays, measure Z. Per-qubit T1 from the same session (the calibration-drift lesson: never trust week-old numbers).
- **Block C (GHZ parity):** prepare (|000⟩+|111⟩)/√2, same delays, 170 kHz virtual detuning on q0, measure all-X. ⟨XXX⟩ envelope = Γ_GHZ(111) on a genuinely entangled carrier.

Envelope fits: damped cosine with the fringe frequency seeded at the virtual detuning (fitted within a bounded window to absorb real qubit-frequency offsets), rate = the envelope exponent. Readout error attenuates amplitudes by (1−2ε)^|S| but is delay-independent, so it does not bias the rates.

## Pre-registered predictions and thresholds

Let σ denote the fit standard error at 4096 shots (simulation calibration: single-pattern rates carry ~2–4·10⁻⁴/µs; the five price readings scatter ~5 % under shot noise alone).

- **P1 (additivity, F70/AT):** each multi-qubit pattern rate equals the sum of its singles within 2σ. Simulation (clean model): excesses ≤ 2.4·10⁻³/µs.
- **P2 (the price, F89d):** the five readings: the three proper pair sums (100|011), (010|101), (001|110), Γ(111) from the product, Γ(111) from the GHZ: coincide. **HOLDS** if relative spread ≤ 10 %; **VIOLATED (correlated dephasing detected)** if ≥ 15 %; in between: inconclusive, repeat with 4× shots before any interpretation.
- **P3 (pattern not carrier, F70):** |Γ_GHZ − Γ_111^product| ≤ 2σ combined. A GHZ excess with clean P2 pairs indicates common-mode noise (superdecoherence), and must then match the common-mode magnitude decoded from P2's fine structure.
- **P4 (T1 accounting, F82):** 1/T2*_j ≥ 1/(2·T1_j) per qubit; Tφ_j = (1/T2*_j − 1/(2T1_j))⁻¹ reported. T1 is local, so the price law holds for the TOTAL rates; no T1 subtraction enters P2.

**Decoding table for a P2 violation** (CORRECTED by the review round; the general formula: the pairing whose singleton is qubit k is elevated by exactly 2·c_ij of the OPPOSITE bond (i,j), so the three pairings form a complete 3-bond covariance tomograph):

| Signature | Reading |
|---|---|
| three proper pair sums EQUAL, each Σγ + 5c; Σ singles at +3c; Γ(111) at +9c (HIGH) | common-mode dephasing of strength c (all-to-all covariance) |
| exactly ONE pairing elevated (+2c), two clean; the elevated one has the singleton k opposite the bond (i,j) | a single correlated bond (i,j) of strength c |
| ONE pairing clean, two elevated | two bonds sharing a qubit (the line's likely case, c₀₁ + c₁₂): the clean pairing's singleton is the SHARED qubit |
| pairings elevated by the boundary rule (middle-qubit pairing highest); Γ(111) AND GHZ clean and LOW | coherent ZZ crosstalk, not dephasing (X_D commutes with Z_iZ_j iff the bond does not cross D's boundary); check against the backend's reported ζ_ij |
| pair sums elevated together with Γ(111) LOW (rates adding in quadrature, not linearly) | LOCAL quasi-static (Gaussian) noise + exponential misfit: an instrument condition, not correlated dephasing; the shape-selection rule must fire before any physics reading |
| pair sums clean, GHZ elevated alone | block-to-block drift or fit artifact in block C; re-fit and check drift before claiming superdecoherence |

**What would falsify what:** a P2 violation does NOT falsify F89d (a theorem); within the Gaussian-dephasing covariance model the five readings coinciding is equivalent to zero pairwise covariance (the tomograph reads each bond separately, so no conspiracy, including anti-correlated bonds, evades it), but OUTSIDE that model the invariance is sensitive to, not equivalent to, locality: coherent ZZ, local non-Markovian shape, and correlated T1 each have their own decoding rows/signatures. Instrument failures (aliasing, fit non-convergence, drift between blocks, shape misfit) are declared as such, not converted into physics claims.

## Simulation record (both stored as JSONs, 2026-07-03)

- Clean heterogeneous model (T2* = 45/88/61 µs, T1 = 180/220/150 µs, readout 2 %): five price readings 0.0490–0.0517/µs around the true P = 0.0500, spread 5.3 % → HOLDS. In-situ T1 recovered 168/223/152 µs.
- Correlated model (common-mode 0.004/µs): spread 32.4 % → VIOLATED fires; additivity excesses carry the |S|² common-mode signature (+0.012/+0.010/+0.008 pairs, +0.034 triple); the proper pairs stay near each other while Γ(111) separates: exactly the decoding table's first row.

## Cost and gating

~33 circuits × 4096 shots ≈ 2 QPU minutes (Batch, ibm_kingston), inside the free monthly tier. Gate order: this document committed → fresh-context adversarial review (two independent lenses, findings folded back into the runner and this spec, recorded below) → fresh calibration pull → Tom's explicit go → the one Batch job.

## Review round (fresh context, before hardware): VERDICTS RECORDED 2026-07-03

Two independent fresh-context reviewers (physics attack; methodology + code). Joint verdict: **NOT ready for hardware as first written; the theory core HELD and got stronger; the analysis pipeline must be rebuilt per the list below.** The deepest shared finding vindicates the review discipline itself: the simulate mode generated data in the fit model's own (wrong) functional form, so the pipeline had validated itself against itself, twice (the fringes AND the T1 asymptote).

**Blockers (all confirmed by derivation):**
- B1: a product state's correlator is a PRODUCT of cosines, not a cosine at the sum frequency (the pair {0,1} carries equal lines at 80 and 20 kHz; the triple four lines at ¼ amplitude). The single-cosine fit biases the four multi-qubit rates pattern-dependently: a near-guaranteed false VIOLATED. Fix: hierarchical fit (singles first, frequencies/phases frozen, correct product model with 3 free parameters per pattern) plus the all-Y quadrature block; simulator must generate the physical signal.
- B2: the T1 analysis assumed relaxation to the maximally mixed state; a transmon relaxes to |0⟩ (⟨Z⟩ = 1 − 2e^(−t/T1)). Fix: 3-parameter fit with free asymptote, physical simulator.
- B3: purely LOCAL Gaussian (quasi-static) dephasing plus an exponential-only fit fakes ~30 % spread (rates add in quadrature). Fix: pre-registered shape-selection (exp·Gauss envelope, AIC) AND the shape-free pointwise form of the law: ln E_D(t) + ln E_D̄(t) = Σ_j ln E_j(t) at every t, whatever the shape; new decoding row (recorded above).

**Majors:** the decoding table's single-bond row was BACKWARDS (one pairing elevated, two clean, not the reverse; corrected above, with the two-bond row added); ZZ crosstalk needs its own row (boundary rule, Γ(111)+GHZ immune; pull the backend's ζ before the run); fringe windows must be absolute (±5 kHz on singles, frequencies FIXED for multi-patterns) with a Nyquist-densified early grid; the verdict bands (HOLDS ≤ 10 % / inconclusive / VIOLATED ≥ 15 %) and the P1/P3 2σ tests must be implemented in code exactly as registered, with proper per-point sigma; the null spread must be calibrated at the COUNTS level (all correlators from the same synthetic shots, per-qubit asymmetric confusion matrix, ≥ 200 seeds) since the readings are correlated through shared shots; blocks must be interleaved (drift); an Aer end-to-end parity run of the actual circuits is mandatory before submission; NaN/failed fits must yield declared instrument failures, never a physics verdict.

**Held under attack:** the pair-sum bookkeeping formula; the common-mode 5c:9c row; the GHZ pattern-not-carrier rate (exact at all orders including T1: no leaked matrix element can dilute it, and every anti-diagonal contaminant decays at the same price); the conventions (no factor-of-2 errors); the bit ordering, JSON round-trip, and delay-granularity handling.

**Gate state:** the runner rebuild per the mandatory list is IN PROGRESS; the re-calibrated thresholds and the Aer parity record will be appended here; hardware submission stays blocked until then, followed by fresh calibration and Tom's explicit go.
