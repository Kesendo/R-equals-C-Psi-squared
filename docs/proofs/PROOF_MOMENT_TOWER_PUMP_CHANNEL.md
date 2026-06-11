# PROOF: the moment-tower pump channel (the device's own damping reads the girth ladder linearly)

**Status:** Tier 1 derived (one-line identities, universal N; verified machine-exact at N = 2, 3). The hardware protocol of §6 **ran on ibm_kingston the same day** (2026-06-11, q149/q13/q9, no entangling gates): the structural law confirmed (double null, row-exact rung-2 identity, girth 2 read from hardware, site tracking, 0.3-5.7% cross-arm reproducibility). The rate layer's first reading (a pump ≤ Γ violation on q13) was corrected within hours by the prep-conditioned re-analysis: the basis-state preparation contains the |0⟩- and |1⟩-branches, so pump and Γ come from the same circuits, the bound holds everywhere in-situ, the margins read the per-qubit thermal population, and the actual finding is minute-scale T1 telegraphing on two of three qubits. The protocol is self-arbitrating. See [experiments/F120_MOMENT_TOWER_KINGSTON.md](../../experiments/F120_MOMENT_TOWER_KINGSTON.md) (the Correction section).
**Date:** 2026-06-11
**Authors:** Thomas Wicht, Claude (Fable 5)
**Builds on:**
- [PROOF_F87_WINDOWED_MONOMIAL_CONVERSE](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) §4 and [`simulations/f87_girth_dichotomy.py`](../../simulations/f87_girth_dichotomy.py): the moment tower t_j(l) = Tr(Z_l H^j), the girth ℓ, and the dichotomy m\* = 2ℓ + 1 when t_ℓ ≠ 0.
- [PROOF_F113_COEFFICIENT_DERIVATION](PROOF_F113_COEFFICIENT_DERIVATION.md): the polarity-asymmetry closed form this channel contains as its first rung.
- [PROOF_F82_T1_DISSIPATOR_CORRECTION](PROOF_F82_T1_DISSIPATOR_CORRECTION.md) and [PROOF_F84_AMPLITUDE_DAMPING](PROOF_F84_AMPLITUDE_DAMPING.md): the pump direction D[σ⁻](I) = +γZ and the net vacuum rate Δγ that weights the channel.
- F117 in [ANALYTICAL_FORMULAS](../ANALYTICAL_FORMULAS.md): the d-leg moment classes, of which this channel reads exactly the d = 1 face.

## Abstract

The girth ladder turned "is this Hamiltonian's palindrome hard, and at which rung?" into moment algebra: the tower t_j(l) = Tr(Z_l H^j) is silent below the girth ℓ, and the first firing rung sets the hardness moment m\* = 2ℓ + 1 with a sum-of-squares coefficient, hard at every γ > 0. This note makes the tower's deg-1 face a *hardware observable*, and the instrument is nothing we add to the device: it is the device's own amplitude damping.

The engine is one line. In the standard noise model (Z-dephasing + amplitude damping), amplitude damping is the unique non-unital piece, and its pump direction is a pure local Z: D[σ⁻_l](I) = +γZ_l (the same (Z_l, I) matrix entry that [F82](PROOF_F82_T1_DISSIPATOR_CORRECTION.md) identified as the dissipator's entire Π²-antisymmetric content), with D[σ⁺_l](I) = −γZ_l, so a thermal device pumps with the **net vacuum rate** Δγ_l = γ↓_l − γ↑_l of [F84](PROOF_F84_AMPLITUDE_DAMPING.md). At the maximally mixed state everything else in the generator annihilates the state, so for any measured observable A,

  **d/dt ⟨A⟩ |_{ρ = I/d} = (1/d) · Σ_l Δγ_l · Tr(A Z_l),**

exactly. Choose A = H^j and the right side is (1/d)·Σ_l Δγ_l·t_j(l): the slope of the j-th energy moment under nothing but the chip's own damping **is the moment tower, read linearly**, rung by rung. The slope is dephasing-blind, evolution-blind, and closes at detailed balance; the first rung is [F113](PROOF_F113_COEFFICIENT_DERIVATION.md) itself, with the exact constant asymmetry = −4^N·slope⟨H⟩; the curvature one order up is *exactly affine* in the generator and fingerprints X/Y-flavored parasites linearly, the complementary lens to F113's Z-drive reader. The certificate is honestly one-sided: a firing rung proves hardness at m\* = 2ℓ + 1; a silent deg-1 tower proves nothing about softness (the k = 4 witness IIXY+ZXZY is silent at every rung and hard at m\* = 11 through its deg-5 class).

The computational anchor is [`simulations/moment_tower_pump_channel.py`](../../simulations/moment_tower_pump_channel.py); the cockpit reading is `framework/diagnostics/f120_moment_tower.py`.

## §1 The law

Work with L = −i[H,·] + Σ_l γ^φ_l D[Z_l] + Σ_l γ↓_l D[σ⁻_l] + Σ_l γ↑_l D[σ⁺_l], σ⁻ = (X+iY)/2, D[c](ρ) = cρc† − ½{c†c, ρ}, on N qubits (d = 2^N).

**Theorem (pump-slope law).** d/dt Tr(A e^{Lt}(I/d)) at t = 0 equals (1/d)·Σ_l Δγ_l·Tr(A Z_l) with Δγ_l = γ↓_l − γ↑_l, for every observable A and every H.

*Proof.* The slope is Tr(A·L(I))/d. The commutator kills I. Dephasing is unital: D[Z_l](I) = Z_l I Z_l − I = 0. Amplitude damping is the non-unital piece: D[σ⁻](I) = σ⁻σ⁺ − ½{σ⁺σ⁻, I} = |0⟩⟨0| − |1⟩⟨1| = +Z per damped site, and D[σ⁺](I) = −Z likewise. So L(I) = Σ_l Δγ_l Z_l. ∎

With A = H^j: **slope_j = (1/d)·Σ_l Δγ_l·t_j(l)**, the deg-1 tower of the girth ladder, weighted by the per-site net damping. The pump is the same object twice over: F82's Π²-antisymmetric entry (the only place the standard dissipator breaks the palindrome's antisymmetric lens) and the generator of this channel. The palindrome-breaking direction and the tower-reading direction are one operator.

## §2 The three blindnesses

All exact, all one-line consequences of the proof above:

1. **Dephasing-blind.** γ^φ_l never enters the slope (unital). A device with any dephasing profile reads the same tower.
2. **Evolution-blind.** The slope does not contain the generator's Hamiltonian at all; only the *measured polynomial* A = H_p^j enters. Changing the evolution H, even entirely, leaves the slope fixed. (This is what pushes the chip's actual Hamiltonian to the curvature, §5, and it is a feature: the slope is a clean baseline that cannot be contaminated.)
3. **Detailed-balance closure.** At γ↓ = γ↑ the pump vanishes and every rung reads zero: the channel is powered exactly by [F84](PROOF_F84_AMPLITUDE_DAMPING.md)'s temperature-independent vacuum component. A device in perfect thermal equilibrium with its drive line has no moment-tower channel; a cold device (γ↑ ≈ 0, the usual transmon situation) has the maximal one.

## §3 Rung one is F113

Within [F113](PROOF_F113_COEFFICIENT_DERIVATION.md)'s scope (H = Σ_l (ω_l/2)Z_l, amplitude damping γ_T1 and pumping γ_pump), two quantities derived in different worlds coincide up to a dimension factor. F113's static Frobenius polarity asymmetry, ‖M₊‖² − ‖M₋‖² = (4^N/2)·Σ_l ω_l(γ_pump,l − γ_T1,l), and this channel's dynamic slope of the energy, slope⟨H⟩ = (1/d)·Σ_l Δγ_l·(ω_l/2)·d = ½·Σ_l Δγ_l ω_l. Comparing:

  **asymmetry = −4^N · slope⟨H⟩**  (exactly; verified to 0.00e+00).

The bridge is worth pausing on. F113's left side is a property of the generator's matrix in operator space, a Frobenius norm imbalance between Π-eigenspaces, with no time in it anywhere. The right side is the initial rate of change of a measured energy. That a static spectral-coordinate imbalance *is* a measurable pump rate (times −4^N) is the channel's claim to the name "F113-analogue" made literal: **the moment-tower channel at rung 1 is F113**, and rungs j ≥ 2 are its ladder, which is exactly what the [connection hunt](PROOF_PI_FACTORS_AS_R_TIMES_D.md) ordered.

## §4 The certificate, honestly one-sided

Scan j = 1, 2, 3, …: the first rung whose slope fires is the girth ℓ (the tower below the girth is identically zero), and then by the [girth dichotomy](../../simulations/f87_girth_dichotomy.py) m\* = 2ℓ + 1, deg = 1, and the palindrome is **hard at every γ > 0** with the sum-of-squares coefficient P_{2ℓ+1,1} = (2ℓ+1)·C(2ℓ,ℓ)·Σ_l t_ℓ². A firing rung is a complete hardness certificate.

Two honesty fences:

- **Silence is not softness.** The channel reads only the d = 1 face of [F117](../ANALYTICAL_FORMULAS.md)'s class decomposition. Hardness can live entirely in higher-leg classes: the k = 4 witness IIXY+ZXZY has t_j(l) = 0 for every j and every site, yet is hard at m\* = 11 through its deg-5 class (p₁₁ = 86507520·γ⁵). A silent deg-1 tower through j_max certifies only that the deg-1 face is dead through m = 2j_max + 1. (The weighted ceiling family behaves the same way: the off-line witness (1,2,1) fires at m\* = 11 in the deg-3 class with its deg-1 tower silent.)
- **Site sums can cancel.** The raw slope reads Σ_l Δγ_l t_ℓ(l), which can vanish accidentally even when t_ℓ ≠ 0 sitewise. Per-site resolution closes the gap: the device's calibration already supplies distinct Δγ_l profiles (per-qubit T1), and deliberately damping one qubit at a time (a weak partial reset on site l) reads t_ℓ(l) individually. The certificate in §6 uses the resolved form.

## §5 The curvature fingerprint, and the complementarity

One order up, the second derivative at the mixed state is

  d²/dt² ⟨A⟩ |_{I/d} = (1/d)·Σ_l Δγ_l·Tr(A·L(Z_l)),

and because the pump itself is generator-independent and L appears exactly once, the curvature is **exactly affine in the generator**, no small-δ expansion needed. For a chip Hamiltonian H_c = H_p + δV against the programmed measurement polynomial A = H_p^j, the dephasing and damping contributions cancel in the difference and

  Δcurvature = (δ/d)·Σ_l Δγ_l·(−i)·Tr(V·[Z_l, A]):

the parasite V is read **linearly against the commutator probes [Z_l, H_p^j]**. The probe set has a flavor: commutators of Z with X/Y-letter strings. So X/Y-flavored parasites with overlap on the probe set are visible at first order (verified: Y₀, Y₀Z₁, Y₁X₂ read to 1.7e-18), while **Z-flavored parasites are exactly invisible** ([Z_l, Z-string] = 0; verified at 0.00e+00). That blind spot is not a defect, it is a division of labor: Z-drives are exactly what F113's balance channel reads linearly, and X/Y-drives are exactly what F113 cannot see (its structural origin: only [Z, σ⁻] ∝ σ⁻ feeds the asymmetry; [X, σ⁻] and [Y, σ⁻] are Hermitian and stay balanced). **The two channels partition the single-site parasite algebra**: F113's balance reads the Z-flavor, the pump curvature reads the X/Y-flavor, and both read linearly with closed-form coefficients.

## §6 The hardware protocol

Everything measurable here is standard hardware vocabulary; nothing needs process tomography.

1. **Preparation.** I/d is prepared in distribution: sample computational basis states |x⟩ uniformly (or twirl), one per shot group. (Our existing hardware arms have used only pure initial states; the basis-average is the one new protocol element, and it is preparation-free per shot.)
2. **Evolution.** Free evolution under the device for a short ramp of times t ∈ {0, δt, 2δt, …} inside the early window (the curvature, dominated by the damping of the pumped Z's, sets the window; the verifier's finite-time block recovers the law from exactly this fit).
3. **Measurement.** ⟨H_p^j⟩ as a Pauli polynomial: for a local programmed H and fixed j the string count is polynomial, and the j = 1, 2 rungs at small N are a handful of settings per time point.
4. **Readout.** The fitted initial slope per rung j, against the prediction (1/d)Σ_l Δγ_l t_j(l) with Δγ_l from calibration (γ↑ ≈ 0 on a cold device, so Δγ_l ≈ 1/T1_l). Site resolution by per-site partial damping where the summed certificate is ambiguous.

The prediction is parameter-free once T1 calibration is read: the slope values, the rung at which they first fire, the detailed-balance null, and the curvature shifts under a deliberately injected X/Y-drive are all closed-form. The natural first target is the girth-2 witness pair (slope⟨H⟩ = 0 exactly, slope⟨H²⟩ = Σ Δγ_l t₂(l)/d ≠ 0): a one-qubit-observable null next to a firing rung on the same device, the same discrimination pattern our [Confirmations](../../simulations/framework/confirmations.py) entries already use.

## §7 Scope, and what is ours

The slope law itself is elementary open-system algebra (non-unital channels pump; unital ones do not), and the d = 1 moments it reads are classical traces of the *programmed* Hamiltonian. What this note banks is the assembly, none of whose pieces needed to be invented but whose composition was unbuilt: the pump direction is F82's antisymmetric entry, its weight is F84's vacuum rate, its first rung is F113's asymmetry with the exact bridge constant −4^N, its tower is the girth ladder's deg-1 face, its silence is bounded honestly by F117's class decomposition, and its curvature splits the parasite algebra with F113 into complementary linear readers. The hardness rung m\* = 2ℓ + 1, until now a statement about power sums of a recentred generator, is the index of the first energy-moment polynomial that responds to the device's own damping.

Fences: deg-1 face only (higher-leg readers, the d = 3 response chain, are named follow-up, not claimed); the protocol of §6 is derived, not yet run; the curvature fingerprint sees X/Y-flavored parasites only (by design, §5); all statements are for the standard noise model (Z-dephasing + amplitude damping; other exotic non-unital channels would add their own pump directions to L(I) and are out of scope).

The verification anchor is [`simulations/moment_tower_pump_channel.py`](../../simulations/moment_tower_pump_channel.py): the pump directions, the law against dense generators, the three blindnesses, the girth certificates with the honest negative control, the three-way F113 bridge (closed form = polarity decomposition = −4^N·slope), the curvature fingerprint with the Z-invisibility, the finite-time protocol fit, and the detailed-balance closure. Every block raises on failure; the process exits 0 only if the whole ledger holds.
