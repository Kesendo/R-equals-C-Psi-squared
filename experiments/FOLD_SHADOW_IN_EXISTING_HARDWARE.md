# The Fold's Shadow in Existing Hardware

**Tier 2** (hardware data finding + corollary-gate result). Date: 2026-07-03. Status: honest record. **Not** a Confirmations-registry entry: the signature is present but degraded, not cleanly measured.

## What this means

Before spending scarce QPU budget on a new run, we asked a simple question: is the thing we want to measure already sitting in data we have? For the fold/price law (F89d) the answer is *almost*. The law's own shape has left a fingerprint in an orphaned hardware run that no document or registry ever indexed, and a second, cleaner fingerprint in a mirror-partner run. But both are shadows: the first is broken by the device noise floor, the second is a single snapshot. The clean law cannot be read off either. So the pre-registered `price_pair` run is genuinely needed, and this note is the honest record of the shadows so they are not lost again.

## The question (the corollary gate)

The fold/price law (F89d, `F89CrossFoldSimilarityClaim`, already Tier-1-derived): under local Z-dephasing, a Pauli-string pattern D and its bitwise complement D̄ decay with rates summing to the total price, Γ(D) + Γ(D̄) = Σⱼ γⱼ, pattern-independent; shape-free form ln E_D(t) + ln E_D̄(t) = Σⱼ ln E_j(t) pointwise. The `price_pair` experiment (`PRICE_PAIR_HARDWARE_PREDICTION.md`) is its hardware test at N=3. Before running it, we swept all ~96 existing hardware result files (plus the repo `data/ibm_*` stores) for the same signature.

## The verdict: not cleanly in existing data

No existing run measures what the law needs, namely complementary pattern envelopes over an idle (Hamiltonian-off) delay grid on shared qubits, with a closing reading (Σγ, or the k=N channel) on the same footing:

- **Hamiltonian-on tomography** (lebensader, iy_yi, zn_mirror, soft_break, y_parity): J = 1 mixes patterns coherently, and there are only 3 to 4 evolution times. No clean single-pattern rate.
- **Single-coherence idle decays** (palindrome, ramsey, tomography, block_cpsi_saturation N=2): exactly one pattern, so a complement is impossible within one run.
- **Cross-run stitch fails by ~6-7×.** The idle single-site rate Σγ_φ from `chain_gamma0` (q12-15, ≈ 0.00505/µs) and the block-channel rates from `block_cpsi_ladder` (≈ 0.13/µs) live on incommensurable effective-γ scales; the block decay is gate/Trotter/SPAM-dominated, roughly 11× the idle γ_φ. The fold cannot be closed across two runs. The law needs Γ(D), Γ(D̄), and Σγ in one session on one qubit set, which is exactly the `price_pair` design.

(Aside: `chain_gamma0`'s multi-qubit trajectory job did not merely lose to fidelity, it crashed on a transpile error, so only its calibration-derived γ survives.)

## The shadow: block_cpsi_ladder (the orphan run)

One run carries the fold's own representation, and no document, experiment writeup, or registry ever mentions it (a `grep` for `block_cpsi_ladder` across all tracked `.py`/`.md`/`.cs` returns zero hits).

In the block-CΨ ladder a coherence |a⟩⟨b| has disagreement k = popcount(a⊕b); its bitwise complement has N−k. So the Hamming channels k and N−k **are** the fold pair D ↔ D̄, and the price law becomes rate(k) + rate(N−k) = rate(N) = 4γN, i.e. the linear ladder rate(k) = 4γk.

- **Simulation (N=5, γ = 0.005/µs):** the law holds exactly. `rung_fits` gives rate(k) = {0.0200, 0.0600, 0.1000} for k = {1, 3, 5}, all R² = 1.0, `checks.all = true`. rate(k) = 4γk to the digit, and rate(5) = 4γ·5 = 0.100 is the price.
- **Hardware (N=4, ibm_kingston q12-15, idle, n_trotter=1):** the complementary pair k=1 and k=3 (= N−1) is present, but the fold is **broken**. Fitted rate(hd=1) = 0.0565/µs (R² = 0.96), rate(hd=3) = 0.0766/µs (R² = 0.81), so **rate(hd=3)/rate(hd=1) = 1.36**, where the linear law demands 3.0. The reason is visible in the raw channel amplitudes:

| t (µs) | hd=1 | hd=3 |
| --- | --- | --- |
| 0 | 0.0854 | 0.0665 |
| 10 | 0.0459 | 0.0196 |
| 20 | 0.0231 | 0.0026 |
| 30 | 0.0125 | 0.0017 |
| 40 | 0.0061 | 0.0021 |
| 50 | 0.0062 | 0.0013 |

The high-Hamming channel, which should decay *fastest* (at 3×), instead falls into the SPAM/readout noise floor (~0.002) by t = 20 µs, so its fitted rate is a noise-floor artifact, not the true 3×. The fold's *shape* is there (two complementary channels on shared qubits, idle dephasing); the fold's *law* is not confirmed (the fast channel is unresolvable at this γ and readout fidelity). That is precisely the failure mode `price_pair` was rebuilt to defeat: low γ, readout mitigation, verdict bands, instrument-failure withholding.

## The second shadow: zn_mirror's mirror-odd observable

In `zn_mirror_ibm_marrakesh_20260429_102824.json` (N=3, J=1, t=0.8, ibm_marrakesh), the two mirror-partner states differ almost entirely in one correlator: ⟨I,X⟩ flips sign with near-equal magnitude, **state_a = +0.649, state_b = −0.674**. Its antisymmetric part a − b = +1.323 dwarfs every other correlator (next largest: Z,Y at −0.36, Y,Z at +0.35). ⟨I,X⟩ is the single mirror-odd observable of the Z_n reflection on this line, a clean and cheap mirror handle that no registry entry names.

## Conclusion

The fold is peeking through existing hardware, in the complementary Hamming channels of an orphan idle run and in a mirror-odd correlator, exactly where the structure predicts. But neither is a measurement of the law: the first is broken by the noise floor (ratio 1.36 vs 3.0), the second is a single snapshot. The clean, low-γ, readout-mitigated, pattern-resolved test remains genuinely necessary, and the shadows are evidence it would succeed. See `PRICE_PAIR_HARDWARE_PREDICTION.md`.

## Provenance

All values verified at source (2026-07-03), external pipeline `.../ibm_quantum_tomography/results/`:

- `block_cpsi_ladder_hardware_ibm_kingston_N4_20260529_014316.json` (rung_fits, rows/by_hd)
- `block_cpsi_ladder_simulate_N5_n2_20260529_011631.json` (rung_fits, checks)
- `zn_mirror_ibm_marrakesh_20260429_102824.json` (expectations state_a/state_b)
- Σγ context: `chain_gamma0_hardware_20260419_105933.json` (gamma_phi_local_per_us; the trajectory job errored)

The fold's *law* is not confirmed on hardware here (it is degraded), so this is a documented shadow, not a confirmation, and it is deliberately kept out of the Confirmations registry.
