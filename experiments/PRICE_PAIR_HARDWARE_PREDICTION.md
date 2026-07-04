# The Price Pair: the F89d fold as a hardware invariant (pre-registered prediction)

**Status:** PREDICTION, written and committed BEFORE any hardware shot (the pre-registration discipline; the fresh-context review round is recorded below before submission). No hardware data exists for this protocol yet (pipeline inventory 2026-07-03: the pattern-pair coherence-decay measurement is absent from all 96 result files). See [`FOLD_SHADOW_IN_EXISTING_HARDWARE.md`](FOLD_SHADOW_IN_EXISTING_HARDWARE.md) for the detailed existing-data check that confirms this: the fold peeks through only as a degraded shadow (the orphan `block_cpsi_ladder` run, fold ratio 1.36 vs the ideal 3.0; and `zn_mirror`'s mirror-odd ⟨I,X⟩), which is why this clean run is genuinely needed.

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

## Rebuild record and pre-hardware amendments (2026-07-03, committed BEFORE submission)

All mandatory fixes are implemented in the rebuilt runner (hierarchical quadrature analysis, counts-level simulator, per-qubit readout mitigation with two calibration circuits, free-asymptote T1, frozen fringe frequencies, dense early grid, interleaved blocks, instrument-failure paths, verdict bands in code). Validation trail, in order:

- The instrument-failure guard fired correctly on the first defective run (aliased phase unwrap): verdict withheld, no physics claim. Fixed by detrending the phase with the programmed fringe.
- Estimator biases found IN SIMULATION and fixed: the log-envelope Jensen bias (debiased quadrature envelopes) and the asymmetric-readout offset (+4-7 % on singles at ε₁−ε₀ = 2 %; fixed by tensor-product readout mitigation from two added calibration circuits). After both: at 10⁶ shots the singles hit the model exactly and all three c readings sit at ±1·10⁻⁴/µs; spread 0.9 %.
- COHERENCE AMENDMENT to the readings: the pair "pattern rate" is defined as the SUM-BRANCH rate (the aligned coherence, γ_i + γ_j + 2c_ij); the branch mean, the product-triple shared-envelope fit and the sum of singles are c-free by construction and serve as the contract cross-checks, excluded from the spread. The FOUR price readings are the three sum-branch pairings and the GHZ. Verified in simulation: a single bond c₀₁ = 0.002 elevates exactly the pairing with singleton q2 and the GHZ by 2c each (decoding table row 2, as corrected); common-mode elevates GHZ 3× the pairings; ZZ = 5 kHz splits the pairings by the boundary rule while GHZ stays exactly clean (its derived immunity, confirmed).
- P2b ADDED as the primary quantitative verdict: the direct per-bond c readings from the branch splitting (c_ij = (Γ₊ − Γ₋)/4 with propagated σ; locality per bond = |c| ≤ 2σ). The spread-based P2 stays as the pre-registered headline with re-calibrated bands.
- NULL CALIBRATION (200 counts-level seeds, 8192 shots): spread median 6.2 %, p90 10.6 %, p99 15.4 % ⟹ AMENDED BANDS: HOLDS ≤ 12 % (~p95), VIOLATED ≥ 18 % (beyond p99.9), inconclusive between. Shots raised to 8192 (~5 QPU min total at ~80 circuits incl. the two calibration circuits).
- AER PARITY (the actual circuits, noiseless, end to end): fringes exactly at the programmed 30.0/50.0/90.0 kHz, all rates ≈ 0 (10⁻⁵..4·10⁻⁴/µs), all c local, decode path clean. The relative-spread verdict now carries a mean floor (P < 0.005/µs ⟹ N.A.) so parity data cannot masquerade as a verdict.

Expected outcomes, restated before the run: the SIMULATOR (local noise by construction) must show "nothing": HOLDS and all c ≈ 0: the contract. The HARDWARE answers with its own voice: the c_ij (and any decoding-table signature). The center is the contract; the finding is the between.

## HARDWARE RECORD (2026-07-04, ibm_marrakesh, job d949n1tgc6cc73fer8sg)

Path [2, 3, 4] chosen from Tom's fresh calibration CSV (05:06Z; min T2echo 181 µs, max readout 1.5 %); the queued Kingston twin (220 pending jobs) was cancelled unbilled; Marrakesh's queue was empty and the 76-circuit Batch (8192 shots, est. ~4.7 QPU min) ran within minutes. Result JSON: `results/price_pair_ibm_marrakesh_20260704_073922.json`.

**Instrument conditions, declared first (per protocol):** the ±5 kHz fringe guard FIRED on q1 (−7.8 kHz) and q2 (+21.5 kHz): real residual qubit detunings, 33 minutes after calibration ⟹ the P2 spread headline is formally WITHHELD. (The quadrature envelopes and the phase-slope measurement are frequency-immune, so the numbers below stay meaningful; the guard exists to keep the pre-registered headline honest, and it did its job.) Shape: the exp·Gauss fit finds REAL quasi-static (Gaussian) content on all three qubits (K² = 1.0/1.1/4.6·10⁻⁴/µs²), so the mixed-sign P1 "EXCESS" flags and the P3 excess land in decoding row 5: the predicted exponential-misfit instrument condition, not correlated-dephasing claims.

**Readings:** in-situ T1 = 316/300/185 µs; T2* ≈ 276/63/56 µs (q0 exceptionally clean; the T2echo ≫ T2* lesson visible on q1/q2). P2b, the primary c-tomograph: c01 = +0.00051 ± 0.00040 (local), c02 = +0.00003 ± 0.00072 (local), **c12 = +0.00126 ± 0.00057 (2.2σ, flagged: a marginal fast-covariance whisper on bond (3,4)**; look-elsewhere over three bonds makes this weak evidence, to be repeated before any claim).

**The finding (the shape-free pointwise law, B3's strong form):** ln E_GHZ(t) − Σ_j ln E_j(t) grows monotonically to **+1.90 at t = 65 µs**: the collective anti-diagonal coherence OUTLIVES the product of its own singles by a large, growing margin. Under local noise of ANY shape this difference equals −Σ Cov_ij(t) ⟹ the slow (quasi-static) dephasing on this line is strongly ANTI-correlated collectively: Var(φ₀+φ₁+φ₂) ≪ Σ Var(φ_j), a partial decoherence-free-subspace behavior of the GHZ coherence. This is the run's genuine hardware voice, visible only in the shape-free form (the fitted exponents alone would have mislabeled it). Status: Tier-2 reading pending a repeat (one run, one line, one device; a constant SPAM offset is excluded by the growth; the small-envelope log-bias is an order below the signal for t ≤ 45 µs).

**Pre-registered expectation vs outcome:** the simulator showed the contract (HOLDS, c ≈ 0, recorded above); the machine revealed (i) real detuning drift within 33 minutes, (ii) a marginal c₁₂ whisper, (iii) anti-correlated slow noise that lets the maximal-disagreement coherence outlive its parts. The center held; the between spoke.

## RUN 2 pre-registration (2026-07-04, before the shot): the repeat on an independent line

Purpose: harden or kill the anti-correlation finding. HYPOTHESIS, stated first: if the anti-correlated slow dephasing is a DEVICE property (shared control/environment), an independent Marrakesh line far from [2,3,4] shows the same growing positive ln E_GHZ − Σ ln E_j; if it is line-specific physics (a local TLS constellation), the new line shows ≈ 0 growth. Either answer is a finding; c₁₂'s 2.2σ whisper gets its independent look-elsewhere test for free. Amendments before the shot (committed): the fringe guard is two-stage (|offset| ∈ (5, 20] kHz = recorded condition, the phase measurement being self-calibrating; > 20 kHz = unwrap validity broken, run voided), and the pointwise shape-free law is computed in-pipeline and stored in the JSON (run 1 computed it offline). Line: [93, 94, 95] from the same 05:06Z calibration (best min T2echo 192 µs; readout 4.9 % on one qubit, handled by the in-job mitigation). Same 76 circuits × 8192 shots ≈ 4.7 QPU min; July total then ≈ 9.4 of the 10 free minutes.

## RUN 2 RECORD (2026-07-04, ibm_marrakesh, job d949tgevtlqs73fu1v30, line [93, 94, 95])

The amended guard behaved as registered: residual detunings −9.4/−8.4 kHz recorded as conditions, run NOT voided; P2 live this time. In-situ T1 = 245/129/242 µs; T2* ≈ 101/71/218 µs.

**Answer 1: the anti-correlation is a DEVICE property.** The in-pipeline pointwise law reads **+1.688 at t = 65 µs** on this chip-distant line (run 1: +1.90): same sign, same magnitude, same growth. The collectively anti-correlated slow dephasing (the GHZ coherence outliving the product of its singles) is Marrakesh-wide, not a three-qubit accident. The pre-registered either/or resolved to: device.

**Answer 2: the c-tomograph caught a real bond, and BOTH rulers agree.** Branch splitting reads **c₀₂ = +0.00253 ± 0.00049 (5.2σ)** on the OUTER pair (93, 95), with the two nearest-neighbour bonds local. Independently, the complement tomograph shows exactly the pairing with singleton q1 (94, the middle) elevated (0.0371 vs ~0.026), which by the corrected decoding row names the OPPOSITE bond (93, 95): the same bond, from a disjoint analysis path. P2 fires VIOLATED (spread 40.4 %) and the decoding attributes it: a fast correlated-dephasing channel between the two outer qubits (shared TLS or control line being the physical suspects), on top of the device-wide anti-correlated slow component (which pulls the GHZ low, consistent with P3's sign). Run 1's marginal c₁₂ (2.2σ, line [2,3,4]) neither repeats nor contradicts: different line, different bonds.

**Bottom line of the two-run campaign (≈ 9.4 free QPU minutes):** the simulator showed the contract; the machine showed (i) real detuning drift on the half-hour scale, (ii) a device-wide ANTI-correlated slow-dephasing component that partially protects the maximal-disagreement coherence, and (iii) one clear correlated fast bond (93,95) read consistently by two independent rulers. The fold framework's price law did exactly what it was built to do: the center held as the contract, and every deviation decoded into a named, localized property of the device's watching. The locality premise of the F-formulas is measurably violated on Marrakesh in two specific, now-quantified ways; both are noise properties of the device, not of the physics the framework models (and the anti-correlated component is the DFS-flavored one: the between outliving its poles).
