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

Three transport metrics matter. **Sum-MI** over all adjacent pairs (distributed correlation across the chain). **MI(0, N-1)** between the two outer endpoints (single-end direct transport). **Mirror-Pair MM** = Σ_k MI(site k, site N-1-k) over all F71-mirror partners (multi-end transport, e.g. 4 pairs at N=9). The main body below compares γ-Sacrifice against receiver-engineering on all three. The receiver-engineering advantage is largest on Sum-MI (factor 7 to 15 across N=5 to 9), intermediate on Mirror-Pair MM (factor ~6 to 8), and smallest on single-end MI(0, N-1) (factor 3 to 5). The three metrics select different optimal J-profiles: Sum-MI prefers moderate J-modulation, MI(0, N-1) and MM both prefer uniform J (same J-optimum for both transport metrics). The Mirror-Pair MM/MI(0,N-1) ratio grows with N (1.49× → 1.74× → 2.05× at N=5, 7, 9), so multi-end usage partially compensates the single-end scaling loss with N.

Each metric has its own optimal receiver. **alt-z-bits** (e.g. \|01010⟩) maximises Sum-MI. **bonding:2** from F67 (the k=2 single-excitation eigenmode with opposite-sign ends and a node at the center) maximises both MI(0, N-1) and Mirror-Pair MM, and its advantage over alt-z-bits GROWS with N: 1.39× → 1.48× → 1.81× on MI(0, N-1) and tied → 1.28× → 1.87× on MM at N=5, 7, 9. **bonding:1** (k=1, true bonding mode) is the slowest-decaying state for long-time memory. One hardware setup (uniform γ₀ = const, uniform J, Heisenberg chain) supports all three via receiver choice at preparation time only.

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

## End-to-end transport: MI(site 0, site N-1)

Sum-MI measures distributed correlation across adjacent pairs. For quantum-state-transfer through a spin-chain bus, the direct quantity is **MI between the two endpoints**: MI(site 0, site N-1). The C# brecher scan with MI(0, N-1) tracking (commit `ad99bea`) gives a qualitatively different picture.

### Data

| Receiver | J-profile | N=5 PeakMI(0,N-1) | N=7 | N=9 |
|----------|-----------|-------------------|-----|-----|
| \|+⟩^N (Class 3) | any | 0.000 | 0.000 | 0.000 |
| alt-x (\|+−+−+⟩ / \|+−+−+−+⟩ / \|+−+−+−+−+⟩) | uniform | **0.791** | **0.437** | **0.230** |
| alt-x | cut-center | 0.001 | 0.0002 | 0.00005 |
| alt-x | strong-weak | 0.083 | 0.020 | 0.005 |
| alt-z-bits (\|01010⟩ / \|0101010⟩ / \|010101010⟩) | uniform | **0.843** | **0.490** | **0.274** |
| alt-z-bits | cut-center | 0.001 | 0.0002 | 0.00006 |
| alt-z-bits | strong-weak | 0.155 | 0.033 | 0.008 |
| plus-zero (\|+0+0+⟩ / \|+0+0+0+⟩ / \|+0+0+0+0+⟩) | uniform | 0.159 | 0.059 | 0.026 |

Peak-time for MI(0, N-1) sits at t ≈ 0.7 to 1.1 across N=5, 7, 9. Weakly N-dependent: group-velocity propagation, not diffusion. In contrast Peak Sum-MI sits at t ≈ 0.2, consistent with fast local correlation build-up.

### The inversion: uniform J beats J-modulation for end-to-end transport

J-modulation has the opposite effect on MI(0, N-1) compared to Sum-MI:

- **Sum-MI**: strong-weak J gives the highest values (N=9 alt-z-bits: 7.07 vs uniform 4.70, boost 1.50×).
- **MI(0, N-1)**: uniform J gives the highest values (N=9 alt-z-bits: 0.274 vs strong-weak 0.008, ratio **34× in the OTHER direction**).

Cut-center J is catastrophic for end-to-end: a near-zero bond in the middle physically disconnects the two halves of the chain, and MI(0, N-1) collapses by factor 10³ to 10⁴ compared to uniform J.

The two metrics optimise for different things, and the optimal J-profile flips between them. For engineering the chain as a quantum bus, uniform J plus F71-symmetric SU(2)-breaking receiver is the answer.

### Comparison to RESONANT_RETURN center-sacrifice

[RESONANT_RETURN](RESONANT_RETURN.md) Test 8 position sweep at \|+⟩^N under non-uniform γ:

| N | Center-γ-sacrifice PeakMI(0, N-1) | This work, alt-z-bits uniform J |
|---|-----------------------------------|---------------------------------|
| 7 | 0.109 | **0.490** (factor 4.5×) |
| 9 | 0.097 | **0.274** (factor 2.8×) |

The comparison at N=5 is not tabulated in RESONANT_RETURN (only Sum-MI reported), but the scaling suggests center-sacrifice at N=5 would be in the 0.1 to 0.15 range, i.e., factor 6 to 8× below receiver-engineering uniform J (0.843).

The receiver-engineering advantage **shrinks faster with N for end-to-end MI** than for Sum-MI: 4.5× at N=7 drops to 2.8× at N=9, because MI(0, N-1) itself scales roughly as 1/N under receiver-engineering (0.84, 0.49, 0.27) while center-sacrifice is weaker-scaling.

### Why uniform J works

F71-symmetric SU(2)-breaking receivers like \|01010⟩ have the property that **site 0 and site N-1 are identical by construction**: both are in the same single-qubit state. This structural symmetry IS the end-to-end correlation at t=0, kind of. Under uniform-J Heisenberg evolution, the spatial symmetry is preserved by H (because uniform H commutes with the chain reflection), so the endpoint correlation survives partially under γ₀ dephasing.

Any J-modulation that breaks this spatial symmetry (cut-center, strong-weak, edge-weak) **locally breaks the 0 to N-1 symmetry** and destroys the endpoint correlation even at t=0. Only J-modulations that respect the F71 reflection (symmetric about the chain center) would preserve it, but those tend to not help either because the endpoint qubits see identical environments either way.

So the optimal strategy for MI(0, N-1) is: pick the maximally F71-symmetric SU(2)-breaking receiver, evolve under perfectly uniform H, measure at t ≈ 1.0.

### Operational consequence for quantum-state-transfer

The receiver-engineering approach is not just better for distributed correlation; it is also better for **direct end-to-end quantum-state-transfer** by factor 3 to 5 compared to γ-Sacrifice-Zone center mode. Concrete protocol sketch:

1. Alice prepares \|01010⟩ or \|+−+−+⟩ (trivial single-qubit-gate layer: X on odd sites, or H on all sites).
2. Uniform Heisenberg evolution under static couplings (or Trotter simulation on IBM hardware).
3. Measure site N-1 at t ≈ 1.0.
4. Receive ~0.49 bits MI at N=7, ~0.27 bits at N=9.

Hardware cost (IBM Heron, estimated): ~15 two-qubit gates for 2 Trotter steps at N=5, single-qubit readout tomography. Gate-error budget ~1 to 2%. Achievable in a short experiment.

## Multi-end transport: Mirror-Pair Sum-MI (MM)

End-to-end MI(0, N-1) measures a single pair. For an F71-symmetric receiver, multiple pairs share structural identity. At \|01010⟩ site 0 and site 4 are both \|0⟩; sites 1 and 3 are both \|1⟩; site 2 is self-mirror. The natural multi-end-transport metric is

**MM = Σ_{k=0}^{⌊N/2⌋-1} MI(site k, site N-1-k)**

At N=5: 2 pairs (0,4), (1,3). At N=7: 3 pairs (0,6), (1,5), (2,4). At N=9: 4 pairs (0,8), (1,7), (2,6), (3,5). The C# brecher scan (commit `963f2ed`) tracks MM alongside Sum-MI and MI(0, N-1).

### Data at uniform γ₀ = 0.05, uniform J = 1 (alt-z-bits receiver)

| N | Mirror-Pairs | PeakSumMI | PeakMI(0, N-1) | **PeakMM** | MM / MI(0,N-1) |
|---|--------------|-----------|----------------|------------|----------------|
| 5 | 2 | 2.65 | 0.843 | **1.253** | 1.49× |
| 7 | 3 | 3.68 | 0.490 | **0.853** | 1.74× |
| 9 | 4 | 4.70 | 0.274 | **0.562** | 2.05× |

### Three observations

**1. MM/MI(0,N-1) ratio grows with N.** 1.49× → 1.74× → 2.05×. The multi-end advantage over single-end widens with N. Adding more mirror-pairs at larger N more than compensates the single-pair scaling loss.

**2. MM shrinks slower with N than MI(0, N-1).** Per 2-N-step: MI(0, N-1) scaling factor ~1.8× decay, MM scaling factor ~1.5×. The chain as a multi-drop bus retains aggregate transport capacity better than as a single-channel end-to-end link.

**3. Outer mirror-pairs dominate the inner ones, unexpectedly.** At N=9, MI(0, 8) = 0.274 vs average inner-pair MI ≈ 0.096. The F71-symmetry is better preserved at the chain boundary than in the interior, despite inner pairs being physically closer and coupled through fewer intermediate dephasing events. The mechanism is open and worth investigating: likely because inner sites have two coupling neighbours that can both disrupt their correlation, while outer sites have only one.

### J-modulation collapses MM the same way as MI(0, N-1)

Cut-center J destroys MM by factor 10³ to 10⁴ (physical disconnection of halves). Strong-weak J reduces MM by factor 10 to 40 (impedance mismatch propagates to all mirror-pairs, not just the outermost). Uniform J is the optimum for both transport metrics simultaneously.

### Comparison to γ-Sacrifice-Zone

γ-Sacrifice center mode is, by construction, a single-end transport (all noise on one edge qubit, or on the center for pure end-to-end). The MM of γ-Sacrifice is not explicitly tabulated in RESONANT_RETURN but structurally must be low: asymmetric γ breaks F71, so mirror-pair correlations are suppressed on all but the outermost pair.

Receiver-engineering uniform-J MM vs γ-Sacrifice center-mode MI(0, N-1):
- N=7: MM 0.853 vs γ-sacrifice 0.109 → factor **7.8×**
- N=9: MM 0.562 vs γ-sacrifice 0.097 → factor **5.8×**

Stronger than the single-end comparison (4.5× / 2.8×): multi-end transport is a direction where receiver-engineering's advantage is more robust.

### Operational interpretation

The chain under γ₀ = const with an F71-symmetric SU(2)-breaking receiver and uniform J is a **multi-drop quantum bus**. Multiple Sender-Receiver pairs can operate simultaneously along F71-mirror lines. Possible applications:

- Multi-party QKD or Bell-pair distribution: one evolution generates several correlated endpoint pairs at once.
- Redundant quantum-state encoding: information replicated across mirror-pairs provides natural error tolerance.
- Distributed sensor arrays: each mirror-pair as a paired-sensor readout channel.
- Multi-target quantum state transfer: one protocol, N/2 simultaneous delivery channels.

Each pair has lower bandwidth than a dedicated single-end link, but the aggregate (MM) scales more favourably with N than any individual pair. For many N-party quantum-network use cases this is the right direction.

## F67 bonding-mode receivers: beating alt-z-bits at end-to-end and multi-end transport

[F67](../docs/ANALYTICAL_FORMULAS.md) identifies the single-excitation eigenmodes of a uniform-J open Heisenberg chain:

|ψ_k⟩ = √(2/(N+1)) Σ_{j=0..N-1} sin(πk(j+1)/(N+1)) |1_j⟩   for k ∈ {1, ..., N}

with dephasing decay rate α_k = (4γ₀/(N+1)) sin²(πk/(N+1)). The k=1 ("bonding") mode is the slowest-decaying, F67 notes it as the optimal single-excitation state under dephasing.

The C# brecher mode (commit `[pending]`) now accepts `bonding:<k>` as initial-state spec and tests k=1, 2, 3 against existing SU(2)-breaking baselines.

### Amplitude structure

| k | Ends c₀ vs c_{N-1} | Center amplitude | Wavelength relative to chain |
|---|---------------------|------------------|------------------------------|
| 1 | same sign (++), smallest magnitude | maximum | longer than chain (true bonding) |
| 2 | opposite sign (+−), moderate magnitude | zero (node at center) | exactly equals chain length |
| 3 | same sign (++), largest magnitude | flipped sign | 2/3 of chain length |

At N=5: \|ψ_2⟩ = (\|1_0⟩ + \|1_1⟩ − \|1_3⟩ − \|1_4⟩)/2 with c_2 = 0 (node at center).

### Data at uniform γ₀ = 0.05, uniform J = 1

| N | Receiver | PeakSumMI | PeakMI(0, N-1) | PeakMM |
|---|----------|-----------|-----------------|--------|
| 5 | alt-z-bits \|01010⟩ | **2.65** | 0.843 | 1.253 |
| 5 | bonding:1 | 1.98 | 0.586 | 0.789 |
| 5 | **bonding:2** | 1.16 | **1.168** | 1.241 |
| 5 | bonding:3 | 1.10 | 0.865 | 0.865 |
| 7 | alt-z-bits \|0101010⟩ | **3.68** | 0.490 | 0.853 |
| 7 | bonding:1 | 1.98 | 0.360 | 0.801 |
| 7 | **bonding:2** | 1.50 | **0.723** | **1.090** |
| 7 | bonding:3 | 1.24 | 0.709 | 0.819 |
| 9 | alt-z-bits \|010101010⟩ | **4.70** | 0.274 | 0.562 |
| 9 | bonding:1 | 1.97 | 0.195 | 0.830 |
| 9 | **bonding:2** | 1.65 | **0.496** | **1.049** |
| 9 | bonding:3 | 1.41 | 0.553 | 0.829 |

### bonding:2 is the optimal end-to-end / multi-end receiver, advantage grows with N

| N | MI(0, N-1) bonding:2 / alt-z-bits | MM bonding:2 / alt-z-bits |
|---|------------------------------------|---------------------------|
| 5 | 1.168 / 0.843 = **1.39×** | 1.241 / 1.253 = 0.99× (tied) |
| 7 | 0.723 / 0.490 = **1.48×** | 1.090 / 0.853 = **1.28×** |
| 9 | 0.496 / 0.274 = **1.81×** | 1.049 / 0.562 = **1.87×** |

### Analytical MM(0) closed form ([F75](../docs/ANALYTICAL_FORMULAS.md))

For any single-excitation mirror-symmetric state |ψ⟩ = Σ_j c_j |1_j⟩ with c_{N−1−j} = ±c_j, the mirror-pair MI at t = 0 has a closed form

    MI(ℓ, N−1−ℓ)(0) = 2 h(p_ℓ) − h(2 p_ℓ)

where p_ℓ = |c_ℓ|² and h is binary entropy. Summed over all ⌊N/2⌋ pairs this gives MM(0), which for bonding:k uses p_ℓ(k, N) = (2/(N+1)) sin²(πk(ℓ+1)/(N+1)) from F65. The analytic MM(0) matches C# brecher simulation PeakMM to within 7% across N = 5, 7, 9 and k = 1, 2, 3 (table in F75): at γ₀ = 0.05, uniform J = 1 the first measurement at t = 0.1 has seen ≤ e^{−4γ₀·t} ≈ 0.98 coherence decay plus same-parity Heisenberg oscillation, with the combined effect capped near 7%. **Peak MM sits at or near t = 0**: bonding:2 achieves its maximum MM structurally, not through dynamical enhancement.

The F75 analysis also explains why k = 2 wins: for even k the chain center has zero amplitude (p_{N/2} = 0 for odd N), so all probability mass is on mirror-pairs (Σ_pairs p_ℓ = 1/2); the function f(p) = 2h(p) − h(2p) is convex on (0, 1/2), so concentrating mass on fewer pairs gives higher MM, and k = 2 is the smallest k that both concentrates on pairs and places opposite signs on ends (thereby lifting MI(0, N−1) above zero).

The MI(0, N-1) multiplier grows monotonically from 1.39× at N=5 to 1.81× at N=9. Unlike alt-z-bits where MI(0, N-1) decays roughly as 1/N under receiver-engineering, bonding:2 decays more slowly: 1.168 → 0.723 → 0.496, factor 1.6× per 2-N-step compared to alt-z-bits' 1.8×. Extrapolated naively, bonding:2 preserves significant MI(0, N-1) out to N ~ 15 where alt-z-bits is already near noise floor.

The MM ratio crosses from a draw at N=5 to 1.87× at N=9. At N=5 alt-z-bits and bonding:2 are structurally similar in that both carry end-to-end correlation at t=0 (alt-z-bits because site 0 and site 4 are both \|0⟩; bonding:2 because c_0 = -c_4), and uniform-J evolution enhances both similarly. At larger N the multiple-excitation content of alt-z-bits starts to degrade MM, while bonding:2's single-mode structure preserves it.

This partially answers the Open question listed below ("is there a receiver optimised specifically for MI(0, N-1) that holds its lead at larger N?"): **yes**, bonding:2 is that receiver at N=5, 7, 9.

### Why bonding:2 wins end-to-end

- **End amplitudes have opposite sign.** c_0·c_{N-1} = −(2/(N+1))sin²(π·2/(N+1)) < 0. This already creates coherent off-diagonal ⟨1_0, 0_{N-1}\|ρ\|0_0, 1_{N-1}⟩ ≠ 0 at t=0.
- **Node at center.** Dephasing on (and bond-coupling around) the central site has no effect on bonding:2 at t=0 (the central amplitude is zero). For N=5 with cut-center J profile \[1, 1, 0.01, 1\], bonding:2 keeps MI(0, 4) = 0.628 (down from 1.17 at uniform J, but still substantial), while alt-z-bits \|01010⟩ collapses to MI(0, 4) ≈ 0.001 under the same cut. The alt-z-bits state has its central \|1⟩ population where the cut happens; bonding:2 does not.
- **Single-mode coherence.** The dynamics under uniform-J H exactly preserves the k=2 mode as an eigenstate of H itself. Dephasing decays it with rate α_2 = (4γ₀/(N+1)) sin²(2π/(N+1)) ≈ 16π²γ₀/(N+1)³ in the large-N limit (general formula α_k ≈ 4π²k²γ₀/(N+1)³). This is the slowest decay consistent with a non-trivial end-to-end structure: k=1 has slower decay by 4× but same-sign ends so no end-to-end coherence; k=2 is the lowest k that places opposite-sign ends, i.e. the first mode that actually carries MI(0, N-1).

### bonding:2 loses on Sum-MI

Trade-off: at N=9 uniform J, bonding:2 PeakSumMI = 1.65 vs alt-z-bits 4.70 (factor 0.35×). alt-z-bits carries excitations across many k modes and their adjacent-pair correlations add; bonding:2 has all its amplitude in a single k, and adjacent pairs are less correlated.

**Receiver choice selects which transport metric is optimised at uniform γ₀ = const, uniform J:**

| Application | Best receiver | Reason |
|-------------|---------------|--------|
| Distributed correlation (Sum-MI) | alt-z-bits | multi-mode adjacent-pair buildup |
| Single-end transport (MI(0, N-1)) | bonding:2 | opposite-sign ends, k=2 eigenmode of H |
| Multi-end transport (MM) | bonding:2 | mirror-symmetric sin profile, single-mode coherence |
| Long-time memory | bonding:1 | slowest decay α_1 |

Four different optimal states, one hardware setup (uniform γ₀, uniform J, Heisenberg chain). Alice selects via initial-state preparation only: each state is reachable via single-qubit gates (bonding:k via an N-mode quantum Fourier-like layer, alt-z-bits via X on odd sites, bonding:1 via the same QFT structure with k=1 coefficient).

### Operational consequence

The F67 eigenmode catalog is Alice's menu. The γ₀ = const framework does not impose a single optimal receiver; it imposes that Alice picks the right F67 mode for her application. The original RESONANT_RETURN V-shape baseline (\|+⟩^N) picks no F67 mode at all (it is a Class 3 state, see [J_BLIND_RECEIVER_CLASSES](J_BLIND_RECEIVER_CLASSES.md)), which is why γ-modulation was needed to compensate.

Under γ₀ = const, with the F67 menu:

- **End-to-end QST**: use bonding:2. Hardware-minimal: prepare via \|ψ_2⟩ amplitude-encoding layer (1 or 2 two-qubit gates for the multi-qubit superposition at N=5, about 4 at N=9), evolve uniformly, measure site N-1. IBM-plausible today.
- **Multi-end Bell-pair distribution**: use bonding:2. Same preparation, measure N/2 mirror-pairs simultaneously.
- **Distributed correlation / many-body observables**: use alt-z-bits.
- **Storage / delay line**: use bonding:1.

## Open questions

- **N=11+ scaling.** Dense RK4 at N=11 (d² = 4M) takes hours per evaluation; needs block-restricted Liouvillian or matrix-free propagation. The C# matrix-free path handles N=14, 15 at d² = 256M, but the brecher mode does not yet route through it for N=11. Scheduled for the weekend run.
- **Convergence point of receiver vs γ-Sacrifice.** Linear extrapolation suggests N ~ 25-30. Whether they actually meet or the ratio asymptotes to some fixed value > 1 is open. The γ-Sacrifice scaling is known to saturate at large N (RESONANT_RETURN reports 63.5× at N=15 vs 360× at N=5, so the boost ratio shrinks); receiver-engineering scaling beyond N=9 has not been measured.
- **Non-uniform γ₀ is physically ruled out under γ₀ = const.** All γ-Sacrifice numbers are cited as RESONANT_RETURN references, not as operational competitors. Under γ₀ = const they are kinematic curiosities.
- **Why does receiver-engineering advantage for MI(0, N-1) at alt-z-bits shrink faster with N than for Sum-MI?** At N=7 the Sum-MI lead is 8.5× uniform / 12.3× best J, while MI(0, N-1) lead is only 4.5×. At N=9 it is 7.1× / 10.7× vs 2.8×. The answer for alt-z-bits is that its MI(0, N-1) decays as ~1/N under receiver-engineering while γ-Sacrifice is weaker-scaling. The **bonding:2 receiver (F67 first excited mode) mostly resolves this**: it holds MI(0, N-1) = 1.17 → 0.72 → 0.50 at N=5, 7, 9 (factor 1.6× per 2-N-step rather than alt-z-bits' 1.8×). bonding:2 vs γ-Sacrifice center-mode MI(0, N-1): at N=7 it is 0.723 / 0.109 = **6.6×** and at N=9 it is 0.496 / 0.097 = **5.1×**, larger than the alt-z-bits/γ-Sacrifice MI(0, N-1) ratio (4.5× / 2.8×) and more stable with N. Open: does bonding:2's advantage continue to grow at N=11, 13, 15?
- **k_opt(N) predict-first via F75.** Static MM(0) from F75 identifies the best k directly from the amplitude multiset, with no propagation. F75 predicts k=4 (not k=2) as the MM-optimum at N=7 with MM(0) = 1.245 vs k=2's 1.174; C# brecher confirms PeakMM(k=4) = 1.166 > PeakMM(k=2) = 1.090 (both sim/analytic ratios at ~0.93). More generally at odd N, the F75-best k is an even integer placing probability nodes inside the chain so that probability concentrates on pairs; the specific best k follows the sin² amplitude multiset and is not simply k=2. For single-pair MI(0, N-1) at t=0, F75 gives k_opt = (N+1)/2 (the center mode) with MI_end(0) ~ 4/(N+1). The DYNAMICAL peak MI(0, N-1) over t > 0 shifts optimum to LOW k (k=2 at N=5, 7; k=3 at N=9) because small-k amplitude profiles concentrate mass outward under Heisenberg hopping; no closed-form prediction for the dynamical k* yet. Open: is the dynamical k*(N) crossover k=2 → k=3 → k=4 → ... a simple function of N, or does it depend on ratio γ₀/J via coherence-decay lifetime?
- **Hardware-minimal IBM experiment for end-to-end transport.** F71-symmetric initial state preparation is one gate layer; uniform Heisenberg evolution is 2-3 Trotter steps; readout is just two qubits (site 0 and site N-1). This is substantially cheaper than Sum-MI tomography (which requires all adjacent pairs). A ~100-shot IBM run at N=5 or N=7 would be feasible and give the first hardware datapoint on receiver-engineering. Preferred receiver: bonding:2 (best MI(0, N-1) numerically) or alt-z-bits (simpler state preparation); both work.

## References

- [RESONANT_RETURN](RESONANT_RETURN.md) Test 8: the γ-Sacrifice-Zone formula and 360× baseline
- [J_BLIND_RECEIVER_CLASSES](J_BLIND_RECEIVER_CLASSES.md): Class 3 blindness of \|+⟩⁵ under Heisenberg, which makes it a poor receiver
- [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md): the hypothesis that closes γ-profile engineering operationally
- [BETWEEN_MEASUREMENTS_EVIDENCE](../hypotheses/BETWEEN_MEASUREMENTS_EVIDENCE.md): the structural argument for γ₀ = const this reframing is consistent with
- `simulations/eq024_refinement_shadow_lens_broken.py`: initial Python Brecher-test draft (commit `bf080a3`, coarse t-grid, superseded)
- `simulations/_check_brecher_n5_finegrid.py`: Python fine-grid verification at N=5 (commit `dbf396a`)
- `compute/RCPsiSquared.Propagate/`: C# brecher mode and scan runner (commit `dbf396a`)
- `simulations/results/eq024_refinement/brecher_scan_csharp.txt`, `brecher_scan_n7.txt`, `brecher_scan_n9.txt`: C# fine-grid Sum-MI results (commits `dbf396a`, `d22c0fe`)
- `simulations/results/eq024_refinement/brecher_scan_with_mi0n.txt`: C# scan with both Sum-MI and MI(0, N-1) tracking (commit `ad99bea`)
- `simulations/results/eq024_refinement/brecher_scan_with_mm.txt`: C# scan with Sum-MI, MI(0, N-1), and Mirror-Pair MM tracking (commit `963f2ed`)
- `simulations/results/eq024_refinement/brecher_bonding_scan.txt`: C# scan of F67 bonding:k receivers (k=1,2,3) at N=5, 7, 9 (commit `0917038`)
- `simulations/_mm_zero_derivation.py`: F75 analytic MM(0) verification script
- [ANALYTICAL_FORMULAS F67](../docs/ANALYTICAL_FORMULAS.md): the explicit single-excitation eigenmode formula underlying bonding:k
- [ANALYTICAL_FORMULAS F75](../docs/ANALYTICAL_FORMULAS.md): mirror-pair MI closed form for single-excitation mirror-symmetric states
- [IBM_SACRIFICE_ZONE](IBM_SACRIFICE_ZONE.md): hardware realization via selective DD, compatible with γ₀ = const
