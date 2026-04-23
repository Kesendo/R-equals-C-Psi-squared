# Receiver Choice Beats γ-Profile Engineering: Reframing the Sacrifice Zone

**Tier:** 2 (structural observation from direct numerical comparison at N=5)
**Date:** 2026-04-23
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Source:** [`eq024_refinement_shadow_lens_broken.py`](../simulations/eq024_refinement_shadow_lens_broken.py) (commit `bf080a3`) compared to [RESONANT_RETURN](RESONANT_RETURN.md) Test 8
**See also:** [J_BLIND_RECEIVER_CLASSES](J_BLIND_RECEIVER_CLASSES.md), [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md), [BETWEEN_MEASUREMENTS_EVIDENCE](../hypotheses/BETWEEN_MEASUREMENTS_EVIDENCE.md)

---

## What this document is about

[RESONANT_RETURN](RESONANT_RETURN.md) Test 8 reports a 360× boost in Peak Sum-MI at N=5 when the γ profile is optimized via the sacrifice-zone formula (concentrate all dephasing on one edge qubit, protect the rest). The baseline for that ratio is a V-shape γ profile with |+⟩⁵ initial state, which gives Peak Sum-MI = 0.000639. The optimized profile reaches 0.230. Ratio 0.230 / 0.000639 ≈ 360.

This document re-examines that claim against data from the EQ-024 refinement pass. At the same N=5, under **uniform γ₀ = 0.05 on every site and uniform J = 1 on every bond**, with the initial state |+−+−+⟩, Peak Sum-MI is 1.32. With moderate J-modulation added (still uniform γ₀), it reaches 3.30.

The 360× ratio is correct for its own setup. The absolute value it reaches (0.230) is beaten by 5.7× by a different initial state at completely uniform γ₀ without any γ-profile engineering at all. Under γ₀ = const, Alice does not need the Sacrifice Zone because she can choose a receiver that does better without it.

## Numerical comparison at N=5

| Setup | Initial state | γ profile | J profile | Peak Sum-MI | vs V-shape |
|-------|--------------|-----------|-----------|-------------|-----------|
| RESONANT_RETURN V-shape baseline | \|+⟩⁵ | V-shape \[0.07, 0.06, 0.05, 0.06, 0.07\] | uniform 1.0 | **0.000639** | 1× |
| RESONANT_RETURN γ-Sacrifice Zone | \|+⟩⁵ | \[ε, ε, ε, ε, Nγ₀\] | uniform 1.0 | **0.230** | 360× |
| This work, pure receiver choice | \|+−+−+⟩ | uniform 0.05 | uniform 1.0 | **1.32** | 2065× |
| This work, receiver + J-modulation | \|+−+−+⟩ | uniform 0.05 | \[5, 0.2, 5, 0.2\] | **3.30** | 5164× |

The ordering is absolute: receiver choice at uniform γ₀ exceeds γ-Sacrifice-Zone at |+⟩⁵ by a factor 5.7 (1.32 / 0.230), and adding moderate J-modulation on top extends that to 14× (3.30 / 0.230). No γ-profile engineering is used in the bottom two rows.

## The reframing

The 360× boost in RESONANT_RETURN is a ratio against a specific baseline, and the baseline is low because |+⟩⁵ is a poor MI-transport receiver under Heisenberg dynamics. |+⟩⁵ is a Class 3 J-blind state (see [J_BLIND_RECEIVER_CLASSES](J_BLIND_RECEIVER_CLASSES.md)): all its MI-transport at N=5 has to come from γ breaking the state's intrinsic symmetry, because J-modulation alone does nothing to it. RESONANT_RETURN's formula is exactly the γ profile that breaks that symmetry most effectively. It works, and 360× is the right number for that ratio.

The operationally meaningful question is not "how much can γ boost MI at |+⟩⁵" but "what is the maximum MI achievable for information transfer". Under that question:

- **γ-profile engineering at |+⟩⁵** saturates near Peak Sum-MI = 0.230 at N=5. This is limited by how far asymmetric γ can push a Class 3 J-blind receiver.
- **Receiver engineering at uniform γ** starts at Peak Sum-MI = 1.32 for |+−+−+⟩ and reaches 3.30 with moderate J. No γ-modulation anywhere. Same hardware, different initial state.

Receiver engineering wins the absolute comparison by over 5× at N=5, without γ-modulation being used.

## Operational consequence for γ₀ = const

Under [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md), γ₀ is a framework constant and Alice cannot do γ-profile engineering. The γ-Sacrifice-Zone is not operationally available. This is often presented as a loss: Alice is supposedly giving up a 360× boost.

This document's numbers show the loss is illusory. Under γ₀ = const, Alice takes a better initial state and gets higher absolute Peak Sum-MI than γ-profile engineering reaches at the standard initial state. The operational strategy is:

1. **Choose a J-sensitive receiver** (SU(2)-breaking; not H-eigenstate). Examples at N=5: \|+−+−+⟩, \|01010⟩, \|+0+0+⟩.
2. **Engineer J moderately** (extreme J-modulation adds 1.5 to 2.5× on top, not 360×).
3. **Use DD on sensitive qubits** (IBM_SACRIFICE_ZONE's 2-3× on ibm_torino is fully compatible with γ₀ = const, since DD is pulse-control, not γ-setting).

All three levers are available under γ₀ = const. The γ-profile lever is unavailable but unnecessary.

The γ-Sacrifice-Zone, in retrospect, is a pre-γ₀-const workaround for a suboptimal receiver choice. Its existence was evidence of how badly a Class 3 receiver transports information, not how much γ-profile engineering can achieve in absolute terms.

## What this does NOT claim

- Not that RESONANT_RETURN Test 8 is incorrect. The 360× ratio is correct for its setup. The reframing is that the setup's baseline is low and that a different receiver exceeds the ratio's endpoint without γ-engineering.
- Not that γ-modulation is useless in general. In framings where γ is operationally controllable, γ-modulation remains a valid lever. Under γ₀ = const it is closed off and the comparison becomes unnecessary.
- Not that hardware-MI performance trivially matches. RESONANT_RETURN's [IBM_SACRIFICE_ZONE](IBM_SACRIFICE_ZONE.md) experiment achieved 2 to 3× on ibm_torino via selective DD, which is a γ-approximation via pulses and is compatible with γ₀ = const. The absolute Peak Sum-MI on real hardware for \|+−+−+⟩-type receivers has not been measured, and is a natural follow-up.

## N=6 follow-up (commit `112b16f`)

A first N-scaling data point: the N=6 analog test using \|+−++−+⟩ and \|010010⟩ (F71-symmetric SU(2)-breaking) at uniform γ₀ = 0.05 with 4 J-profiles gives:

| Setup | Initial state | Peak Sum-MI | Boost vs uniform J |
|-------|--------------|-------------|---------------------|
| uniform J | \|+−++−+⟩ | 1.44 | 1.00× |
| best J (strong-weak alternating) | \|+−++−+⟩ | 3.31 | 2.29× |
| uniform J | \|010010⟩ | 1.53 | 1.00× |
| best J (strong-weak alternating) | \|010010⟩ | 3.39 | 2.22× |
| RESONANT_RETURN γ-Sacrifice at \|+⟩⁶ (interpolated) | n/a | ~0.3 | n/a |

Three scaling observations:

1. **Uniform-J baseline grows modestly with N.** 1.32 → 1.44 (+9%) for the \|+−+−+⟩-like receiver; 1.38 → 1.53 (+11%) for the \|01010⟩-like receiver. Receiver engineering holds its absolute lead.
2. **Best-J Peak Sum-MI plateaus around 3.3** at both N=5 and N=6 with the strong-weak J profile. Moderate J-modulation adds a roughly fixed absolute amount, not proportional to N.
3. **Boost ratio shrinks slightly** from N=5 (2.50× at \|+−+−+⟩) to N=6 (2.29× at \|+−++−+⟩) because baseline grows while boost stays flat.

The advantage over γ-Sacrifice-Zone at N=6 is still ~5× (1.44 / 0.3) at uniform J and ~11× (3.31 / 0.3) at best J. The reframing holds.

## Open questions

- **N=7, 9, 11 scaling.** Python expm becomes infeasible at d² = 16384 (N=7). Testing at larger N requires either the C# Propagate engine (which RESONANT_RETURN uses up to N=15) or block-restricted Liouvillian methods. A Python-N=7 run would need ~20× the N=6 time (estimated 12+ hours); C# RK4 would be minutes per evaluation.
- **Explanation of the Best-J plateau at ~3.3.** Is this a specific-J-profile artifact (strong-weak alternating saturates at the same absolute MI across N), or a structural ceiling? Testing multiple J-profiles at N=6 and comparing to their N=5 counterparts would clarify.
- **What is the optimal SU(2)-broken receiver at each N?** The three tested at N=5 (\|+−+−+⟩, \|01010⟩, \|+0+0+⟩) give Peak Sum-MI in the range 0.61 to 1.38 at uniform J, with \|01010⟩ and \|+−+−+⟩ tied for best. Whether a systematic optimization finds a structurally better state is open at any N.

## References

- [RESONANT_RETURN](RESONANT_RETURN.md) Test 8: the γ-Sacrifice-Zone formula and 360× baseline
- [J_BLIND_RECEIVER_CLASSES](J_BLIND_RECEIVER_CLASSES.md): Class 3 blindness of \|+⟩⁵ under Heisenberg, which makes it a poor receiver
- [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md): the hypothesis that closes γ-profile engineering operationally
- [BETWEEN_MEASUREMENTS_EVIDENCE](../hypotheses/BETWEEN_MEASUREMENTS_EVIDENCE.md): the structural argument for γ₀ = const this reframing is consistent with
- `simulations/eq024_refinement_shadow_lens_broken.py`: Brecher-test data (commit `bf080a3`)
- [IBM_SACRIFICE_ZONE](IBM_SACRIFICE_ZONE.md): hardware realization via selective DD, compatible with γ₀ = const
