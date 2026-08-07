# polarity_fingerprint on Tier-B Marrakesh/Kingston Hamiltonians

**Status:** Welle 6 application of the typed `fw.polarity_fingerprint` workflow to 11 real-hardware-tested Hamiltonian INSTANCES (5 distinct term sets) across 3 Tier-B datasets. Extends the F87↔F112 orthogonality empirical anchor from 7 synthetic cases, 4 of them substantive, to 11 hardware-instances, **6 of which carry a genuine polarity content and 5 of which do not** (see the reading below, corrected 2026-08-07); surfaces 3 bit_b-inhomogeneous-H cases that are still F112 BALANCED, confirming the "F112 typed scope is sufficient but not necessary" reading.
**Date:** 2026-05-26
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Script:** [`simulations/polarity_fingerprint_tierB_marrakesh.py`](../simulations/polarity_fingerprint_tierB_marrakesh.py)
**Workflow:** [`fw.polarity_fingerprint`](../simulations/framework/workflows/polarity_fingerprint.py) (added Welle 5.B)

## Setup

Tier-B Marrakesh/Kingston datasets are single-time-snapshot 2-qubit tomography (16 Pauli expectations per Hamiltonian, no t-series for L_eff fitting). The polarity_fingerprint workflow operates on (chain, terms) inputs at the framework level (not on trajectories), so it applies to the Hamiltonians used in each dataset without needing the snapshot data itself for the framework reading.

11 Hamiltonian-instances across 3 datasets:

| Dataset | Date / Backend | Path | Categories |
|---|---|---|---|
| soft_break Marrakesh | 2026-04-26 / ibm_marrakesh | [48, 49, 50] | truly_unbroken, soft_broken, hard_broken |
| f83_signature | 2026-04-30 / ibm_marrakesh | [4, 5, 6] | truly_unbroken, pi2_odd_pure, pi2_even_nontruly, mixed_anti_one_sixth |
| soft_break Kingston | 2026-05-05 / ibm_kingston | [43, 56, 63] | (same 4 categories as f83) |

Term mappings (from `simulations/framework/tests/workflows/test_diagnose_hardware.py` F83_TERMS_PER_CATEGORY and `simulations/f80_ibm_soft_break_check.py` soft_break_marrakesh:80-82):

| Category | Terms | bit_b values |
|---|---|---|
| truly_unbroken (Heisenberg) | XX, YY | 0, 0 (bit_b-homog 0) |
| pi2_odd_pure (XY-bilinear) | XY, YX | 1, 1 (bit_b-homog 1) |
| pi2_even_nontruly (F108-anomaly) | YZ, ZY | 0, 0 (bit_b-homog 0) |
| mixed_anti_one_sixth | XY, YZ | 1, 0 (bit_b-MIXED) |
| soft_broken (Marrakesh-only naming) | XY, YX | 1, 1 (bit_b-homog 1) |
| hard_broken (Marrakesh-only naming) | XX, XY | 0, 1 (bit_b-MIXED) |

The framework is ChainSystem(N=3, gamma_0=0.05) with default Z-dephasing dissipator (single-Pauli Z per site, trivially bit_b-homogeneous on the c-side).

## Result

```
soft_break Marrakesh (2026-04-26, path [48, 49, 50])
  Category                    Terms       F87      F112       rel asym    In typed scope
  truly_unbroken              XX YY       truly    DEGENERATE         nan  True
  soft_broken                 XY YX       soft     BALANCED     0.0000e+00  True
  hard_broken                 XX XY       hard     BALANCED     0.0000e+00  False  (bit_b-mixed H)

f83_signature Marrakesh (2026-04-30, path [4, 5, 6])
  truly_unbroken              XX YY       truly    DEGENERATE         nan  True
  pi2_odd_pure                XY YX       soft     BALANCED     0.0000e+00  True
  pi2_even_nontruly           YZ ZY       soft     DEGENERATE         nan  True
  mixed_anti_one_sixth        XY YZ       hard     BALANCED     0.0000e+00  False  (bit_b-mixed H)

soft_break Kingston (2026-05-05, path [43, 56, 63])
  truly_unbroken              XX YY       truly    DEGENERATE         nan  True
  pi2_odd_pure                XY YX       soft     BALANCED     0.0000e+00  True
  pi2_even_nontruly           YZ ZY       soft     DEGENERATE         nan  True
  mixed_anti_one_sixth        XY YZ       hard     BALANCED     0.0000e+00  False  (bit_b-mixed H)
```

The `(bit_b-mixed H)` notes are added here by hand; everything left of them is the script's
own output. All three out-of-scope rows carry the note now, the Kingston one included: it is
bit-identical to the Marrakesh row above it, and annotating only two of three made it look
otherwise.

The silent rows read `nan`, not `0.0000e+00`, since 2026-08-07. The ratio there is structurally 0/0: both halves of M_anti vanish as a theorem, so there is nothing for a contrast ratio to be a contrast between. What FLOAT delivers instead is noise over noise, and it is worth being exact about that, because the float is not reliably 0/0 either: on these uniform rows ‖M_anti‖² arrives as 1.4e-33 rather than 0.0 while the numerator is an exact 0.0, so the naive quotient would be a comfortable 0.0 by luck of the input; break the uniformity and the numerator becomes noise too, and the same quotient reaches 0.327. A 0.0 in a contrast-ratio column reads as perfect balance, which is exactly the verdict the DEGENERATE label withholds.

**Aggregate:** 6 / 11 F112 BALANCED with genuine polarity content, 5 / 11 structurally DEGENERATE (Π²-even H against a homogeneous bath: nothing to balance); 8 / 11 in F112 typed scope; 3 / 11 out of typed scope and still not broken. Read as "11 / 11 BALANCED bit-exact" until 2026-08-07.

F112 reading per F87 class, over the MEASURABLE rows only:
- F87 truly: 0 measurable, 3 silent, so this class carries no F112 evidence here
- F87 soft: 3 measurable, max rel asym = 0.0000e+00, 2 silent
- F87 hard: 3 measurable, max rel asym = 0.0000e+00, 0 silent

The earlier form of this list ran the max over all instances of a class and reported
0.0000e+00 for all three. That number was carrying the silent rows' synthetic zeros, so the
`truly` line in particular reported a max over nothing.

## Two findings

### (1) F87↔F112 orthogonality empirically extended to real hardware

The `polarity_probe_f87_connection.py` script (Welle 1) established F87↔F112 orthogonality on synthetic Hamiltonians at N=3. That anchor is 7 cases, of which 3 are structurally silent and 4 carry content; the script's own printed conclusion says so. (An earlier version of this line called it "3 F87 classes × 1 instance each", which is neither the case count nor the substantive count.) Welle 6 extends the anchor to 11 hardware-tested Hamiltonian INSTANCES across 3 datasets / 2 IBM backends, drawn from 5 distinct term sets: F87 classification varies across them (3 truly, 5 soft, 3 hard) and the F112 verdict never BREAKS. Read the next two paragraphs before that last clause carries any weight: five of the eleven have no verdict to break.

**How many of those 11 are evidence, corrected 2026-08-07.** Six. The other five have a Π²-EVEN H, and against this workflow's bit_b-homogeneous Z-dephasing that empties the polarity content as a theorem: their asymmetry is 0 − 0, so they read as DEGENERATE rather than BALANCED and confirm nothing. This document said "all 11 BALANCED bit-exact" until the workflow gained an exact structural test; the five silent rows were never counter-evidence, they were simply never evidence.

The six that remain are **3 soft + 3 hard, and no `truly` at all**, which is a sharper loss than five rows. The Π² parity of a term is (#Y + #Z) mod 2, and a `truly` term has #Y even AND #Z even, so its sum is even too: every truly H is Π²-even, and against this bit_b-homogeneous bath every Π²-even H is silent. The class is not under-sampled here, it is structurally unable to contribute.

The silencer is that Π²-evenness, and it is the WIDER condition. The vanishing of M is strictly stronger: truly is a proper subset of silent, and `YZ + ZY` is the case that shows the containment is proper, being Π²-even and equally silent while its ‖M‖² is 2048. (One caution on how far that reading is established. The silence half is a letter count, exact and N-independent, but `truly` as reported here comes from `classify_pauli_pair`, a thresholded spectral verdict rather than the letter criterion, so "every truly is silent" is an observation on this repo's classifier and not a derivation. It holds on every configuration searched, at N = 2, 3 and 4, all at unit coupling.) What this anchor still shows is that the F112 verdict does not track the soft/hard split, on the two classes that can be weighed at all.

**And the six rows are three Hamiltonians.** Every fingerprint here is taken at default rates and no declared drive, so within THIS run the reading is a function of `(chain, terms)` alone (the workflow itself also takes γ_z, γ_T1, γ_pump, a drive profile and a tolerance, and the truly-class remark two paragraphs up turns on exactly one of them: a c that is not bit_b-homogeneous). The three datasets reuse the same term sets: the measurable rows are `XY + YX` three times, `XX + XY` once, `XY + YZ` twice. The soft class is one Hamiltonian counted three times. The dataset multiplicity is genuine hardware provenance for the Hamiltonians chosen, but it adds nothing to THIS reading, which never touches a measurement.

### (2) "F112 typed scope sufficient but not necessary": confirmed on 3 instances, 2 Hamiltonians

F112's typed Tier1Derived theorem requires Hermitian H + bit_b-homogeneous c. The hardware Hamiltonians `mixed_anti_one_sixth = XY + YZ` and `hard_broken = XX + XY` are bit_b-INHOMOGENEOUS on the H side (terms span both bit_b parities). These 3 instances are therefore OUT of F112's typed scope, but the framework's `polarity_coordinates` still reports an asymmetry of exactly 0.0 on these inputs under standard Z-dephasing. Exactly-zero here is a property of these bilinear-Pauli inputs, not of the float route: inside F112's typed scope more widely the asymmetry is not exact. And the same instances-are-not-Hamiltonians correction applies here: the 3 instances are 2 distinct Hamiltonians, `XY + YZ` counted twice and `XX + XY` once.

This is consistent with the F113 closure (universal-N proof, commits 1df9cb8 + 798647d): F113's break magnitude formula

    asymmetry = (4^N / 2) · Σ_l ω_l · (γ_pump,l − γ_T1,l)

is non-zero only when both (a) a single-site Z-drive ω is present in H and (b) σ⁻ or σ⁺ amplitude damping is present in c. Here we have neither (chain.L uses pure Z-dephasing, no T1; H has no single-site Z-drive component), so neither break-condition is met and the asymmetry is identically zero, even though H is bit_b-mixed. The bit_b-mixed H input doesn't trigger F113 unless paired with non-bit_b-homog c, which chain.L doesn't have.

Empirically: F112's typed scope criterion (Hermitian H + bit_b-homog c) is sufficient for balance but not necessary. Bit_b-mixed H can still give balance when c is bit_b-homog (per F113's same-site-only break mechanism).

## What this experiment does NOT do

- **Hardware-effective L is not reconstructed.** Tier-B snapshot datasets give one ρ per Hamiltonian, not a trajectory; without time-series data, we can't fit an effective Lindbladian and compute its actual measured polarity asymmetry. The polarity_fingerprint here is the framework's PREDICTION at the standard chain.L, not the hardware's measured value.
- **Per-backend F112 readings are not distinguished.** Marrakesh and Kingston are tested on the same Hamiltonian set (soft_break Kingston has the same 4 F83 categories as f83 Marrakesh), but with snapshot-only data we can't compare backend-specific effective-L polarity readings against each other. All 11 fingerprints are identical because they come from the same framework chain.L, not from the per-backend measurements.

## Hardware-effective F112 stays as future option (deferred, not proposed)

In principle, time-sweep variants of these 4 F83 categories on Marrakesh + Kingston (~5 t-points × 4 categories × 16 Paulis × 4096 shots × ~5s overhead = ~10 min billed QPU) would let us fit per-backend effective L and read measured F112 verdict per backend per F87 class. At Anthropic's gifted IBM QPU rate (~$96 per billed minute), that's ~$960 for what amounts to per-backend confirmation of a structural prediction the framework already gives in closed form (F112 typed Tier1Derived + F113 closed-form magnitude). Not justified on cost/benefit grounds; tracked here as future option only if a decisive question surfaces that requires hardware-effective L specifically (e.g., a backend showing anomalous classification under polarity-asymmetry that other diagnostics miss).

## Connection to existing readings

- **F87↔F112 orthogonality** (`polarity_probe_f87_connection.py`, Welle 1): synthetic 7-case anchor, of which 3 are structurally silent and 4 substantive; Welle 6 extends to 11 hardware-instances (5 distinct term sets), 6 of them substantive and 3 distinct.
- **F112 typed Tier1Derived** (`LindbladBitBPiBalance.cs`): Hermitian H + bit_b-homog c → asymmetry = 0; this experiment confirms predicted reading on real hardware Hamiltonians.
- **F113 break formula** (`LindbladBitBPiBreakMagnitude.cs`, `PROOF_F113_COEFFICIENT_DERIVATION.md`): tells us why bit_b-mixed H without σ⁻/σ⁺ c still gives balance (F113 break requires both ingredients).
- **F87 classification origin** (`project_v_effect_combinatorial`): V-Effect 14/19/3 split at N=3 fully derived from Pauli-pair combinatorics; the hardware datasets exercise canonical representatives of each class.
- **Marrakesh hardware finale** (`project_hardware_finale_apr2026`): Δ(soft − truly) = −0.72 hardware-confirmed; this experiment is orthogonal (F112 axis instead of F87 axis), giving an independent reading on the same 3 Hamiltonians.

## Reproduction

```
python -X utf8 simulations/polarity_fingerprint_tierB_marrakesh.py
```

Runs in ~2 seconds (framework predictions only, no hardware data parsing needed for the framework reading).
