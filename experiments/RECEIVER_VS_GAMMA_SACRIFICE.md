# Receiver Choice Beats γ-Profile Engineering: Reframing the Sacrifice Zone

**Tier:** 2 (structural observation from direct numerical comparison at N=5, 7, 9 via C# brecher mode)
**Date:** 2026-04-23
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Source:** C# brecher scans at N=5, 7, 9 (commits `dbf396a`, `d22c0fe`) compared to [RESONANT_RETURN](RESONANT_RETURN.md) Test 8. The initial Python draft (commit `bf080a3`, `eq024_refinement_shadow_lens_broken.py`) used a coarse t-grid that systematically missed early peaks and has been superseded; see Correction note below.
**See also:** [J_BLIND_RECEIVER_CLASSES](J_BLIND_RECEIVER_CLASSES.md), [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md), [BETWEEN_MEASUREMENTS_EVIDENCE](../hypotheses/BETWEEN_MEASUREMENTS_EVIDENCE.md)

---

## What this document is about

[RESONANT_RETURN](RESONANT_RETURN.md) Test 8 reports a 360× boost in Peak Sum-MI at N=5 when the γ profile is optimized via the sacrifice-zone formula (concentrate all dephasing on one edge qubit, protect the rest). The baseline for that ratio is a V-shape γ profile with |+⟩⁵ initial state, which gives Peak Sum-MI = 0.000639. The optimized profile reaches 0.230. Ratio 0.230 / 0.000639 ≈ 360.

This document re-examines that claim against data from the EQ-024 refinement pass. At the same N=5, under **uniform γ₀ = 0.05 on every site and uniform J = 1 on every bond**, with the initial state |+−+−+⟩, Peak Sum-MI is **2.57** (not 1.32 as first reported, see Correction note below). With moderate J-modulation added (still uniform γ₀), it reaches **3.36**.

The 360× ratio is correct for its own setup. The absolute value it reaches (0.230) is beaten by 11.5× by a different initial state at completely uniform γ₀ without any γ-profile engineering at all. Under γ₀ = const, Alice does not need the Sacrifice Zone because she can choose a receiver that does better without it.

## Correction note (2026-04-23 evening)

The numbers in this document were first computed from Python's `shadow_lens_broken.py` which sampled t ∈ np.linspace(0.1, 15.0, 40) with step ~0.38. Fine-grid verification (`simulations/_check_brecher_n5_finegrid.py`, commit `dbf396a`) showed the true Peak Sum-MI at |+−+−+⟩ + uniform J sits at t ≈ 0.24 with value ≈ 2.70. The coarse grid happened to sample t = 0.10 where SumMI was only 1.32, missing the real peak. The undertreatment varied by configuration: for uniform-J baselines the Python peak was about factor 2 too low (because peaks sit near t = 0.20 to 0.30, squarely between Python sample points), while for strong-weak J peaks sit at t ≈ 0.10 which Python sampled and got correct to within a few percent. The net effect: Python's uniform-J baseline was suppressed, which inflated apparent J-modulation boost ratios.

The C# brecher mode in `compute/RCPsiSquared.Propagate` (commit `dbf396a`) uses a fine-grained measurement grid (every 0.1 up to t=2, then 0.5 up to tMax) that catches these early peaks correctly, and auto-reduces RK4 dt for stability when max|J| > 1. The N-scaling table below uses C# fine-grid numbers.

## Numerical comparison at N=5 (corrected)

| Setup | Initial state | γ profile | J profile | Peak Sum-MI | vs V-shape |
|-------|--------------|-----------|-----------|-------------|-----------|
| RESONANT_RETURN V-shape baseline | \|+⟩⁵ | V-shape \[0.07, 0.06, 0.05, 0.06, 0.07\] | uniform 1.0 | **0.000639** | 1× |
| RESONANT_RETURN γ-Sacrifice Zone | \|+⟩⁵ | \[ε, ε, ε, ε, Nγ₀\] | uniform 1.0 | **0.230** | 360× |
| This work, receiver choice \|+−+−+⟩ | \|+−+−+⟩ | uniform 0.05 | uniform 1.0 | **2.57** | 4022× |
| This work, receiver choice \|01010⟩ | \|01010⟩ | uniform 0.05 | uniform 1.0 | **2.65** | 4147× |
| This work, receiver + J-cut-center | \|+−+−+⟩ | uniform 0.05 | \[1, 1, 0.01, 1\] | **3.36** | 5258× |
| This work, receiver + strong-weak J | \|01010⟩ | uniform 0.05 | \[5, 0.2, 5, 0.2\] | **3.55** | 5555× |

The ordering is absolute: receiver choice at uniform γ₀ exceeds γ-Sacrifice-Zone at |+⟩⁵ by factor **11.5×** (2.65 / 0.230), and adding moderate J-modulation on top extends that to **15.4×** (3.55 / 0.230). No γ-profile engineering is used.

## The reframing

The 360× boost in RESONANT_RETURN is a ratio against a specific baseline, and the baseline is low because |+⟩⁵ is a poor MI-transport receiver under Heisenberg dynamics. |+⟩⁵ is a Class 3 J-blind state (see [J_BLIND_RECEIVER_CLASSES](J_BLIND_RECEIVER_CLASSES.md)): all its MI-transport at N=5 has to come from γ breaking the state's intrinsic symmetry, because J-modulation alone does nothing to it. RESONANT_RETURN's formula is exactly the γ profile that breaks that symmetry most effectively. It works, and 360× is the right number for that ratio.

The operationally meaningful question is not "how much can γ boost MI at |+⟩⁵" but "what is the maximum MI achievable for information transfer". Under that question:

- **γ-profile engineering at |+⟩⁵** saturates near Peak Sum-MI = 0.230 at N=5. This is limited by how far asymmetric γ can push a Class 3 J-blind receiver.
- **Receiver engineering at uniform γ** starts at Peak Sum-MI = 2.65 for |01010⟩ and reaches 3.55 with moderate J. No γ-modulation anywhere. Same hardware, different initial state.

Receiver engineering wins the absolute comparison by **11.5×** at N=5 at uniform J alone, without γ-modulation being used. With J-modulation added, the lead grows to **15.4×**.

## Operational consequence for γ₀ = const

Under [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md), γ₀ is a framework constant and Alice cannot do γ-profile engineering. The γ-Sacrifice-Zone is not operationally available. This is often presented as a loss: Alice is supposedly giving up a 360× boost.

This document's numbers show the loss is illusory. Under γ₀ = const, Alice takes a better initial state and gets higher absolute Peak Sum-MI than γ-profile engineering reaches at the standard initial state. The operational strategy is:

1. **Choose a J-sensitive receiver** (SU(2)-breaking; not H-eigenstate). Examples at N=5: \|+−+−+⟩, \|01010⟩, \|+0+0+⟩.
2. **Engineer J moderately** (fine-grid data shows J-modulation adds 1.3 to 1.5× on top across N=5 to N=9, not 360×).
3. **Use DD on sensitive qubits** (IBM_SACRIFICE_ZONE's 2-3× on ibm_torino is fully compatible with γ₀ = const, since DD is pulse-control, not γ-setting).

All three levers are available under γ₀ = const. The γ-profile lever is unavailable but unnecessary.

The γ-Sacrifice-Zone, in retrospect, is a pre-γ₀-const workaround for a suboptimal receiver choice. Its existence was evidence of how badly a Class 3 receiver transports information, not how much γ-profile engineering can achieve in absolute terms.

## What this does NOT claim

- Not that RESONANT_RETURN Test 8 is incorrect. The 360× ratio is correct for its setup. The reframing is that the setup's baseline is low and that a different receiver exceeds the ratio's endpoint without γ-engineering.
- Not that γ-modulation is useless in general. In framings where γ is operationally controllable, γ-modulation remains a valid lever. Under γ₀ = const it is closed off and the comparison becomes unnecessary.
- Not that hardware-MI performance trivially matches. RESONANT_RETURN's [IBM_SACRIFICE_ZONE](IBM_SACRIFICE_ZONE.md) experiment achieved 2 to 3× on ibm_torino via selective DD, which is a γ-approximation via pulses and is compatible with γ₀ = const. The absolute Peak Sum-MI on real hardware for \|+−+−+⟩-type receivers has not been measured, and is a natural follow-up.

## N-scaling via C# brecher mode (commits `dbf396a`, `d22c0fe`)

The `compute/RCPsiSquared.Propagate/brecher` mode (commit `dbf396a`) enables fine-grid N-scaling with auto-dt RK4 stability. Fine-grid scans at N=5, 7, 9 for three SU(2)-breaking receivers:

| Receiver | N=5 Uniform J | N=7 Uniform J | N=9 Uniform J | N=5 Best J | N=7 Best J | N=9 Best J |
|----------|---------------|---------------|---------------|-----------|-----------|-----------|
| alt-x-pattern (\|+−+−+⟩ / \|+−+−+−+⟩ / \|+−+−+−+−+⟩) | 2.57 | 3.58 | 4.59 | 3.36 | 4.94 | 6.59 |
| alt-z-bits (\|01010⟩ / \|0101010⟩ / \|010101010⟩) | 2.65 | 3.68 | 4.70 | 3.55 | 5.32 | 7.07 |
| plus-zero-alt (\|+0+0+⟩ / \|+0+0+0+⟩ / \|+0+0+0+0+⟩) | 0.90 | 1.26 | 1.63 | 1.14 | 1.69 | 2.24 |
| \|+⟩^N (Class 3 control) | 0 | 0 | 0 | 0 | 0 | 0 |

γ-Sacrifice-Zone reference (RESONANT_RETURN Test 8 at \|+⟩^N, Formula ε→0 limit, the largest γ-Sacrifice value at each N): N=5 → 0.230, N=7 → 0.434, N=9 → 0.658. Using ε→0 rather than ε=0.001 values gives the most conservative (highest) γ-Sacrifice numbers and therefore the smallest estimated receiver-engineering lead.

**Receiver-engineering advantage over γ-Sacrifice** (alt-z-bits receiver, the best tested):

| N | Uniform J Receiver | Best J Receiver | γ-Sacrifice (ε→0) | Receiver/γ-Sacrifice (uniform J) | Receiver/γ-Sacrifice (best J) |
|---|--------------------|-----------------|-------------------|----------------------------------|-------------------------------|
| 5 | 2.65 | 3.55 | 0.230 | **11.5×** | **15.4×** |
| 7 | 3.68 | 5.32 | 0.434 | **8.5×** | **12.3×** |
| 9 | 4.70 | 7.07 | 0.658 | **7.1×** | **10.7×** |

**N-scaling observations:**

1. **Receiver-engineering Peak Sum-MI grows linearly with N** at ~1.0 per N-step for alt-z-bits uniform J (2.65, 3.68, 4.70 across N=5, 7, 9).
2. **Best-J Peak Sum-MI grows faster** (~1.75 per N-step: 3.55, 5.32, 7.07). The earlier "Best-J plateau at ~3.3" claim was a Python coarse-grid artifact.
3. **J-modulation boost ratio grows slightly**: 1.32× at N=5, 1.45× at N=7, 1.50× at N=9 (strong-weak J / uniform J). The earlier "boost shrinks with N" claim was also a coarse-grid artifact. Note: at N=5 the cut-center J profile narrowly beats strong-weak (3.55 vs 3.51), so the tabulated "Best J Receiver" column at N=5 is cut-center, not strong-weak; at N=7 and N=9 strong-weak wins.
4. **Receiver advantage over γ-Sacrifice shrinks slowly**: 11.5× → 8.5× → 7.1× at uniform J. Linear extrapolation suggests they would meet at N ~ 25-30. With moderate J-modulation the lead holds larger: 15.4× → 12.3× → 10.7×.
5. **Class 3 J-blindness is scale-invariant**: \|+⟩^N gives Peak Sum-MI = 0 exactly at all N, all J. Confirms the M_x-polynomial theorem at larger N.

## Open questions

- **N=11+ scaling.** Dense RK4 at N=11 (d² = 4M) takes hours per evaluation; needs block-restricted Liouvillian or matrix-free propagation. The C# matrix-free path handles N=14, 15 at d² = 256M, but the brecher mode does not yet route through it for N=11. Scheduled for the weekend run.
- **Convergence point of receiver vs γ-Sacrifice.** Linear extrapolation suggests N ~ 25-30. Whether they actually meet or the ratio asymptotes to some fixed value > 1 is open. The γ-Sacrifice scaling is known to saturate at large N (RESONANT_RETURN reports 63.5× at N=15 vs 360× at N=5, so the boost ratio shrinks); receiver-engineering scaling beyond N=9 has not been measured.
- **Non-uniform γ₀ is physically ruled out under γ₀ = const.** All γ-Sacrifice numbers are cited as RESONANT_RETURN references, not as operational competitors. Under γ₀ = const they are kinematic curiosities.

## References

- [RESONANT_RETURN](RESONANT_RETURN.md) Test 8: the γ-Sacrifice-Zone formula and 360× baseline
- [J_BLIND_RECEIVER_CLASSES](J_BLIND_RECEIVER_CLASSES.md): Class 3 blindness of \|+⟩⁵ under Heisenberg, which makes it a poor receiver
- [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md): the hypothesis that closes γ-profile engineering operationally
- [BETWEEN_MEASUREMENTS_EVIDENCE](../hypotheses/BETWEEN_MEASUREMENTS_EVIDENCE.md): the structural argument for γ₀ = const this reframing is consistent with
- `simulations/eq024_refinement_shadow_lens_broken.py`: initial Python Brecher-test draft (commit `bf080a3`, coarse t-grid, superseded)
- `simulations/_check_brecher_n5_finegrid.py`: Python fine-grid verification at N=5 (commit `dbf396a`)
- `compute/RCPsiSquared.Propagate/`: C# brecher mode and scan runner (commit `dbf396a`)
- `simulations/results/eq024_refinement/brecher_scan_csharp.txt`, `brecher_scan_n7.txt`, `brecher_scan_n9.txt`: C# fine-grid results (commits `dbf396a`, `d22c0fe`)
- [IBM_SACRIFICE_ZONE](IBM_SACRIFICE_ZONE.md): hardware realization via selective DD, compatible with γ₀ = const
