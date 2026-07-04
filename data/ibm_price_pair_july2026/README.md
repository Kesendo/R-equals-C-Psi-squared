# The price-pair campaign, ibm_marrakesh, 2026-07-03/04

The complete data of the four-run pre-registered campaign testing the F89d fold at the ρ level
(the price law Γ(D) + Γ(D̄) = Σγ as a device invariant, i.e. the LOCALITY premise under the
F-formulas). Pre-registration, review round, all run records and the campaign verdict:
[experiments/PRICE_PAIR_HARDWARE_PREDICTION.md](../../experiments/PRICE_PAIR_HARDWARE_PREDICTION.md).
Plain-words arc: [reflections/ON_THE_STOLEN_MEMORY.md](../../reflections/ON_THE_STOLEN_MEMORY.md).
Runner (external pipeline): `AIEvolution.UI/experiments/ibm_quantum_tomography/run_price_pair.py`
(modes `--simulate/--null/--aer/--hardware/--zztest/--analyze`; the hardware JSONs here re-analyze
with `--analyze FILE`).

| File | What |
|------|------|
| `price_pair_ibm_marrakesh_20260704_073922.json` | RUN 1, line [2,3,4], job d949n1tgc6cc73fer8sg: fringe guard fired (detuning drift −7.8/+21.5 kHz), c-tomograph c₁₂ 2.2σ (did not recur), pointwise law +1.90 at 65 µs |
| `price_pair_ibm_marrakesh_20260704_075304.json` | RUN 2, line [93,94,95], job d949tgevtlqs73fu1v30: pointwise law +1.688 (device-wide), c₀₂ = +0.00253 ± 0.00049 read by both rulers (downgraded after run 3: did not recur) |
| `price_pair_ibm_marrakesh_20260704_082618.json` | RUN 3 (W discriminator), job d94abhvu62ks7395p4ig: H1 refuted, W pairs +1.505/+0.946/+1.718 = the 3:2:3 boundary-rule fingerprint of coherent ZZ |
| `price_pair_zztest_ibm_marrakesh_20260704_083938.json` | RUN 4 (conditional Ramsey), job d94ajjtgc6cc73fes7bg: ζ₀₁ = −3.92 ± 0.14 kHz, ζ₁₂ = −3.64 ± 0.10 kHz, ζ₀₂ = +0.02 ± 0.12 kHz: the ZZ attribution confirmed deterministically |
| `price_pair_null_20260703_214747.json` | the 200-seed counts-level null calibration at 8192 shots (median 6.2 %, p90 10.6 %, p99 15.4 %) that set the amended verdict bands (HOLDS ≤ 12 %, VIOLATED ≥ 18 %) |
| `price_pair_aer_20260703_214803.json` | the Aer end-to-end parity record of the actual circuits (noiseless: fringes exactly 30/50/90 kHz, rates ≈ 0) |
| `price_pair_zztest_sim_20260704_083506.json` | sim validation of the ZZ discriminator (5 kHz model read as −4.2/−4.1/−0.0 incl. the declared T1 attenuation ≈ 0.84) |

**Campaign verdict, one line:** the dephasing is local (all covariances null in clean sessions);
the one systematic non-locality is coherent nearest-neighbour ZZ ≈ −3.9/−3.6 kHz (null
next-nearest), first seen as a Gaussian-masked pseudo-anti-correlation, fingerprinted by the W
run, then measured by conditioning; hardware noise model for the F-formulas = local Z-dephasing
(quasi-static + Markovian) + local T1 + coherent NN ZZ ≈ 4 kHz, no correlated bath.
