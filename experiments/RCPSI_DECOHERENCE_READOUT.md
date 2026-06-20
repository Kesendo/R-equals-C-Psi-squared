# R=CΨ² as a decoherence readout: the deviation from Born is a γ-meter

**Status:** Tier 2/3. The leading-order deviation is Tier-1 (F94); the invertibility/monotonicity of the
readout is gate-verified numerically (not yet an analytic theorem).
**Date:** 2026-06-20
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 4.8)
**Grounded in:** [F94](../docs/ANALYTICAL_FORMULAS.md) (`Δ = P_lindblad/P_unitary − 1 = (4/3)J²γt³`,
Tier-1 derived, [`PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md`](../docs/proofs/PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md)).
**Builds on:** [`BORN_RULE_MIRROR.md`](BORN_RULE_MIRROR.md) (the generalized Born rule `R_i = C_i·Ψ_i²`,
Tier 2/3); complementary to [`GAMMA_AS_SIGNAL.md`](GAMMA_AS_SIGNAL.md) (a different readout, below).
**Verifiers:** [`_rcpsi_readout_f94.py`](../simulations/_rcpsi_readout_f94.py) (the F94 substrate, all gates),
[`_rcpsi_readout_rabi.py`](../simulations/_rcpsi_readout_rabi.py) (the textbook Rabi qubit),
[`_rcpsi_as_readout.py`](../simulations/_rcpsi_as_readout.py) (the refuted N=2 |++⟩ null).

## What this is about

A textbook says: measure a qubit and the outcome probabilities follow the Born rule, `P_i = |⟨i|ψ⟩|²`.
Open the system to a dephasing environment and the measured populations drift a few percent off that
Born prediction. The usual name for that drift is *noise*. This note says the opposite: the drift is
*signal*. The amount by which the open-system population deviates from the closed-system Born value is a
clean, monotone, invertible function of the decoherence rate γ, so **measuring the deviation reads γ from
inside** — no separate probe, no environment access. The system's own departure from the textbook rule
*is* the measurement of its environment.

## The mechanism: R = C·Ψ²

Write the generalized Born rule ([`BORN_RULE_MIRROR.md`](BORN_RULE_MIRROR.md)) per outcome:

> `R_i = C_i · Ψ_i²`,  with `Ψ_i² = ⟨i|U ρ₀ U†|i⟩` (closed-system Born population, unitary)
> and `R_i = ⟨i|ρ_open(t)|i⟩` (open-system population, Lindblad), so `C_i = R_i / Ψ_i²`.

`C_i` is the per-outcome "mirror quality": Born is the `C ≡ 1` perfect-mirror limit. The content — what
makes it a readout — is that `C_i` is a **monotone invertible function of γ**. Physically: the Hamiltonian
rotates coherence into population and back; dephasing intercepts the coherence mid-rotation, changing how
much amplitude lands in outcome `i`. So `C_i − 1` measures how much the environment has interfered with the
coherent rotation. The readout only has content where `H` *makes* that coherence↔population transfer; in a
stationary/perfect-mirror state (`C ≡ 1`) there is nothing to read.

## Grounded in F94 (why it inverts)

The invertibility is not assumed; it falls out of a proven F-formula. F94 (Tier-1 derived) is exactly this
deviation, for the dominant outcome `|00⟩` of pair (0,2) of `|0+0+⟩` on the N=4 Heisenberg ring under
Z-dephasing:

> `Δ = C − 1 = P_lindblad/P_unitary − 1 = (4/3)·Q²·K³ = (4/3)·J²·γ·t³`  (deep perturbative regime).

`Δ` is **linear in γ** (at fixed J, t), hence strictly monotone, hence locally invertible, with the
explicit inverse `γ ≈ (C − 1)/((4/3)·J²·t³)`. F94 supplies the proven leading order; the readout is its
operational reading — *the proven Born-deviation, used as a γ-meter.*

## Gate-verified result

On F94's own substrate ([`_rcpsi_readout_f94.py`](../simulations/_rcpsi_readout_f94.py), reusing F94's exact
functions — no convention guessing):

- **The deviation is F94's:** in the deep regime, `Δ/((4/3)J²γt³) → 0.999` as `t → 0` (and `Δ ∝ γ^0.98`).
- **Monotone and broadly invertible:** `C(γ)` is strictly monotone increasing over the whole tested range
  `γ ∈ [0.01, 10]` (`Q = J/γ` from 100 down to 0.1) — **no turnover** at `t = 1`. The readout is not
  confined to a thin perturbative window.
- **Inversion works:** a hidden γ is recovered from the measured `C` (exact `C(γ)` calibration) to `~1e-5`;
  the F94 leading-order inverse tracks it in the deep regime.

The textbook Rabi qubit ([`_rcpsi_readout_rabi.py`](../simulations/_rcpsi_readout_rabi.py); `H = J·X`, `|0⟩`,
Z-dephasing) shows the same: `C_|1⟩(γ) = R/sin²(Jt)` strictly monotone, hidden γ recovered to `~1e-6`.

## Scope and honesty

- **Where it works:** wherever `H` rotates coherence into populations and dephasing intercepts (the Rabi
  qubit; F94's `|0+0+⟩`-ring). **Where it does not:** the stationary / perfect-mirror limit — the N=2 `|++⟩`
  Heisenberg case ([`_rcpsi_as_readout.py`](../simulations/_rcpsi_as_readout.py)) has stationary populations,
  `R = Ψ²` exactly, `C ≡ 1`, nothing to read. That take-1 is a documented null, not a success.
- **Tier.** The leading-order deviation `(4/3)J²γt³` is Tier-1 (F94). Global monotonicity / invertibility is
  **gate-verified numerically** over the tested range, not an analytic theorem; a Zeno turnover at extreme
  dephasing (`γ ≫ J`) is physically expected but does not appear up to `Q = 0.1` here. So the readout is
  Tier 2/3: a proven seed (F94) with a numerically-verified invertibility.
- **A convention caveat learned here.** F94's Hamiltonian is `(J/4)·Σ(XX+YY+ZZ)` per bond, not
  `J·Σ(XX+YY+ZZ)`; a naive build with the wrong normalization produced a *spurious* turnover at `γ ≈ 0.4`.
  Anchoring on F94's own functions removed it. (The deviation value is convention-sensitive; the
  invertibility thesis is not.)

## Relation to GAMMA_AS_SIGNAL (a different readout)

[`GAMMA_AS_SIGNAL.md`](GAMMA_AS_SIGNAL.md) is the flagship "decoherence is signal" result, but a *different*
readout: it recovers the **spatial γ-profile** (the per-site γ vector) from a full observable set
(purities + pair CΨ + MI) via the palindrome-full-rank Jacobian / SVD channel (15.5 bits). This note is the
**scalar dual**: a single γ, read from one outcome's Born deviation `C = R/Ψ²`, tied directly to the
generalized Born rule rather than to the γ-Jacobian. The two are complementary faces of the same principle —
the chain reads its own environment from the inside.

## What this adds

The framing "`C_i = R_i/Ψ_i²` is a monotone invertible decoherence readout — the deviation from Born *is*
the measurement" is new (not previously typed; only `BORN_RULE_MIRROR`'s `R=C·Ψ²` decomposition and F94's
deviation scaling existed). The contribution is the *operational inversion* and its grounding: F94's proven
`Δ ∝ γ` is what makes the namesake formula a γ-meter.
