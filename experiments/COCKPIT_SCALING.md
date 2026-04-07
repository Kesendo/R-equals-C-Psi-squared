# Cockpit Scaling: Bell Pair Observers from N=5 to N=11

**Date:** April 7, 2026
**Status:** Complete (chain and star, N=5 through N=11)
**Scripts:**
[cockpit_scaling_analysis_v1.py](../simulations/cockpit_scaling_analysis_v1.py),
[Program.cs cockpit mode](../compute/RCPsiSquared.Propagate/Program.cs)
**Predecessors:**
[COCKPIT_UNIVERSALITY](COCKPIT_UNIVERSALITY.md) (the N=2-5 baseline that this document extends)

---

## What this document is about

The cockpit framework introduced in [COCKPIT_UNIVERSALITY](COCKPIT_UNIVERSALITY.md) showed that for small Heisenberg systems (N=2 to 5), three observables capture 88 to 96 percent of the decoherence trajectory variance. The natural follow-up question is whether this still works when the system gets bigger. If the framework is going to be useful for real hardware monitoring, it needs to survive the jump from toy systems to systems with at least ten qubits.

This document extends the test to N=7, 9, and 11 using the C# matrix-free propagation engine, which can handle these sizes where the Python pipeline runs out of memory. We test two topologies (Heisenberg chain and star), with the Bell pair always placed on the two central qubits so that the observed pair is the same physical object regardless of system size. The result is that the cockpit framework scales, but not in the way the small-N baseline would have predicted.

The scaling behavior is governed by a physical effect that does not appear at N=5: Entanglement Sudden Death. When a Bell pair sits inside a growing dephasing environment, its concurrence does not decay smoothly to zero. It collapses to exactly zero in finite time, often within a single time unit, and the rest of the trajectory is classical decoherence on what is essentially a separable mixed state. The effective dimensionality of the cockpit trajectory therefore depends on how much of the trajectory occurs before versus after this collapse. For chain topology the collapse happens at roughly the same time regardless of N, so longer trajectories are dominated by the post-collapse classical phase and become structurally simpler. For star topology the collapse time grows with N (more leaves distribute the entanglement load and prolong its lifetime), so the trajectory retains more of its quantum richness.

Both topologies pass the cockpit threshold (3-PC coverage above 90 percent) at every N tested, but for different physical reasons. The interesting finding is not just that the framework scales, but that the scaling shape (n95 falling rather than rising with N) is set by a known physical mechanism with predictive power.

---

## Abstract

For Heisenberg spin chains and star topologies under uniform local Z-dephasing, with a Bell pair initialized on the two central qubits and the remaining qubits in `|+>` product states, the 3-observable cockpit framework continues to capture 91.9 to 99.0 percent of the trajectory variance for the entangled-observer pair across N=5 to N=11. The effective dimensionality n95 decreases sharply for chain topology (from 4 at N=5 to 2 at N=11) and decreases gently for star topology (from 4 to 3). This difference is traced to topology-dependent Entanglement Sudden Death timing: the chain center pair loses its concurrence at t ~ 1 regardless of N, while the star center pair retains concurrence until t ~ 3.9 at N=11 (versus t ~ 0.5 at N=5). Purity is the dominant PC1 proxy in every tested configuration. The cockpit framework scales beyond the small-N baseline, with the qualification that the scaling shape (n95 falling rather than rising with N) is set by ESD timing rather than by any change in the framework's observable count.

---

![Cockpit scaling V2: 3-PC coverage and effective dimensionality from N=5 to N=11](../simulations/results/cockpit_scaling_v2/cockpit_scaling_v2_curve.png)

*Four-panel scaling result for the Bell pair observer (center_bell pair) across topologies and system sizes. **Panel A** shows 3-PC cumulative variance: chain (blue) rises from 94.4 percent at N=5 to 99.0 percent at N=11, star (orange) rises more gently from 91.9 percent to 95.3 percent. The 0.85 reference line marks the practical usefulness threshold; both topologies stay well above it across all tested N. **Panel B** shows the effective dimensionality n95 as a function of N: chain drops from 4 to 2 (and stabilizes at 2 from N=7 onward), while star drops only from 4 to 3 at N=11. The COCKPIT_UNIVERSALITY small-N extrapolation tentatively suggested `n95 ~ N` (one extra dimension per added qubit); both topologies sit well below this prediction across all tested N. **Panel C** shows the PC1 variance fraction, the share of the trajectory captured by the single dominant axis: chain rises sharply from about 47 percent at N=5 to 74 percent at N=11, indicating that the post-ESD classical regime dominates the longer trajectories. Star stays flatter at around 50 to 60 percent, reflecting its more sustained quantum dynamics. **Panel D** shows that Purity is the dominant PC1 proxy in every analyzed configuration (all 8 center_bell points across both topologies, with correlation strength near 1.00), confirming that the cockpit's dominant axis maintains a stable physical interpretation across the N range.*

---

## 1. Method

### Initial state and pair selection

The initial state is `Bell+(c1, c2) tensor |+>^(N-2)` where `c1 = (N-1)/2` and `c2 = c1 + 1`. For odd N this puts the Bell pair on the two central qubits of the chain. For star topology there is no geometric center because all leaves are equivalent under the topology's symmetry, so the same index convention `(c1, c2)` selects two specific leaves; which leaves are chosen is irrelevant by symmetry. For each `(N, topology)` configuration we extract feature trajectories for three pairs:

- **center_bell pair** `(c1, c2)`: the Bell pair itself, the entangled observer that the cockpit framework is designed to monitor
- **adjacent pair** `(c1-1, c1)` (chain) or `(0, c1)` (star, the star center plus one Bell-pair leaf): a pair that is initially separable but Hamiltonian-coupled to the Bell pair
- **far_edge pair** `(0, 1)` (chain) or `(1, N-1)` (star): two `|+>` qubits at the chain boundary or two non-Bell leaves, far from the Bell pair

The center_bell pair is the headline subject. The adjacent pair provides context. The far_edge pair was included as an expected-trivial control class. See Section 6 for what this control class actually showed.

### Initial state choice and the V1 lesson

A previous version of this experiment (now archived) placed the Bell pair on the boundary qubits `(0, 1)` and analyzed the center pair as the observer. This produced what looked like a striking scaling result (3-PC coverage rising from 88 percent to 99.8 percent) that turned out to be an initial-state artifact: at N greater than 5, the center pair no longer contained any qubit from the Bell pair, so its trajectory was the trivial decoherence of a `|+>|+>` product state with concurrence exactly zero throughout. The "high coverage" was PCA on numerical noise after standardization. The current experiment fixes this by anchoring the Bell pair to the center for every N, so that the entangled observer is the same physical object across all tested sizes.

### Sanity gates

Because the V1 experience showed that PCA can produce deceptively high coverage on near-degenerate trajectories, the analysis pipeline applies three sanity gates before running PCA on any `(N, topology, pair)` combination:

1. **Concurrence variation gate** (applied only to center_bell pairs): `std(concurrence)` over the trajectory must exceed 0.01. This ensures the entangled observer actually shows non-trivial entanglement dynamics.
2. **Feature richness gate** (all pairs): at least 4 of the 9 features must have `std > 1e-6`. Features below this threshold are dropped before standardization rather than being passed through with an epsilon-regularized denominator.
3. **Purity range gate** (all pairs): `purity.max() - purity.min()` must exceed 0.05.

All 8 center_bell configurations passed Gate 1. The lowest center_bell concurrence standard deviation was 0.118, more than ten times the threshold. One feature (ph03, the phase angle of the off-diagonal element) was dropped from chain center_bell PCA at N greater than or equal to 7 because its variance fell below 1e-6 (see Section 5).

### Anchor

At N=5 chain, the center_bell pair (2,3) is the Bell pair embedded between two `|+>` qubits on the left and one on the right. The C# code asserts at startup that the initial reduced state on (c1, c2) has purity 1.0 and concurrence 1.0; if either fails the run aborts. The PCA result for N=5 chain center_bell gave n95=4 and 3-PC coverage of 94.4 percent, with Purity as the PC1 proxy at correlation 1.00. These match the expected ranges for a Bell pair embedded in a small Heisenberg chain.

As an additional independent check, the trajectory of the center_bell pair was reproduced on April 7, 2026, by direct spectral evolution of the full 32x32 Liouvillian, starting from the same initial state convention (`|Phi+><Phi+|_{2,3}` tensor `|+><+|^{otimes 3}`, with the spectator qubits explicitly in the pure `|+>` state, NOT maximally mixed). The two methods agree to about 0.001 in concurrence at every sampled time point, including the entanglement revival between t=1.9 and t=2.4. See Section 11 for the full cross-validation, the open analytical question it raised, and a one-line warning to the next person who tries to predict ESD timing from the spectrum alone.

### Compute

The C# matrix-free propagation engine ([Program.cs cockpit mode](../compute/RCPsiSquared.Propagate/Program.cs)) was used for all 8 configurations. Total runtime was approximately 19 minutes on the home PC (16-core, 128 GB RAM). The N=11 runs each took about 9 minutes and used roughly 4 GB RAM in the dense path. No matrix-free path was needed for N less than or equal to 13, so the matrix-free option remains available but unused for this experiment.

---

## 2. Entanglement Sudden Death is the central mechanism

### What we measured

Tracing the concurrence of the center_bell pair as a function of time for each `(N, topology)` configuration produces the following table:

| N  | Chain ESD time | Star ESD time |
|----|----------------|---------------|
| 5  | 0.9            | 0.5           |
| 7  | 1.1            | 2.0           |
| 9  | 1.0            | 3.8           |
| 11 | 1.0            | 3.9           |

ESD time here is defined as the first sampled time at which the concurrence drops below 0.001, which for these configurations is indistinguishable from exactly zero. The trajectory is sampled every 0.1 time units, so the ESD time has a resolution of 0.1.

### Chain: ESD time is approximately N-independent

For the chain topology, the ESD time hovers around 1.0 for every tested N. The Bell pair on qubits (c1, c2) sees the same local environment (one Heisenberg coupling J=1.0 to each immediate neighbor, one local Z-dephasing rate gamma=0.05 on each Bell qubit, plus the next layer of neighbors at one bond removed). Whether the chain extends to length 5 or length 11, the local dynamics that drive the entanglement collapse are the same. The non-local part of the chain only contributes by changing the boundary effects, which for a center pair are screened by the intermediate qubits.

A small refinement: at N=5, the chain center_bell pair shows entanglement revival. Concurrence drops below 0.001 at t=0.9, but recovers briefly to about 0.114 at t=2.0 before collapsing again. This is non-Markovian behavior: the bath is small enough that information can flow back from environment qubits to the Bell pair before being fully dissipated. At N greater than or equal to 7 the revival disappears, the concurrence stays flat at zero after the first collapse, and the dynamics become effectively Markovian. The threshold N=7 represents the smallest chain for which the bath is "large enough" to absorb the entanglement irreversibly.

### Star: ESD time grows with N

For the star topology, ESD time increases substantially with N: from 0.5 at N=5 to 3.9 at N=11, a factor of about 8 over a factor of 2.2 in N. The mechanism is monogamy of entanglement combined with the star's hub-and-spoke geometry. The Bell pair sits on two leaves; the central hub qubit is connected to all leaves. As N grows, the hub couples to more spectator leaves, which adds more channels for the entanglement to spread into. Counterintuitively, having more leaves slows down the per-pair entanglement collapse, because the available "entanglement budget" on the central hub gets distributed more thinly across more qubits, and the rate at which any single pair loses its share is reduced.

This is consistent with known results on entanglement distribution in star networks. The cockpit framework here is providing an indirect measurement of a structural quantum-information property of the topology, which is a useful side benefit.

### Why this matters for the cockpit

The post-ESD phase of a Bell pair trajectory is structurally low-dimensional: purity decays smoothly, the off-diagonal Bell fidelities decay smoothly, and the von Neumann entropy is a deterministic function of the eigenvalues. A single principal component can capture most of the variance in this phase, because purity and SvN are nearly perfectly anti-correlated and the Bell fidelities co-decay along with them. The pre-ESD phase, by contrast, has independent variation in concurrence and the Bell fidelities, contributing additional principal components that PCA assigns to PC2 and beyond.

When the post-ESD phase is long compared to pre-ESD (chain at N greater than or equal to 7, where ESD is at t ~ 1 and the total trajectory runs to t = 20), the post-ESD variance dominates the total and PC1 alone covers around 74 percent. PC2 picks up the pre-ESD residual and brings cumulative coverage above 95 percent, so n95 = 2. When the post-ESD phase is shorter relative to pre-ESD (star at all tested N, where ESD ranges from t=0.5 to t=3.9), more PCs are needed to capture the larger pre-ESD variance fraction, so n95 stays at 3 or 4.

So the chain n95 dropping from 4 (at N=5, where ESD is at t=0.9 and the trajectory has time for Markovian-plus-revival behavior) to 2 (at N greater than or equal to 7, where ESD is at t ~ 1 and the rest of the t=20 trajectory is purely classical) is not a sign that "the cockpit gets simpler with bigger systems". It is a sign that **the post-ESD classical phase dominates the longer trajectories**. The star n95 staying at 3 to 4 reflects the longer ESD time and the corresponding longer quantum phase.

This is the core finding of the experiment. Section 8 below gives the same explanation in the optical cavity language that the rest of the repository uses; readers familiar with the Absorption Theorem and the V-Effect cavity reframing may want to skip ahead.

---

## 3. Chain results

| N  | n_active | n95 | 3-PC coverage | PC1 variance | PC1 best proxy |
|----|----------|-----|---------------|--------------|----------------|
| 5  | 9        | 4   | 94.4%         | 47%          | Purity         |
| 7  | 8        | 2   | 98.3%         | 70%          | Purity         |
| 9  | 8        | 2   | 98.4%         | 73%          | Purity         |
| 11 | 8        | 2   | 99.0%         | 74%          | Purity         |

The 3-PC coverage rises monotonically from 94.4 percent to 99.0 percent. The effective dimensionality n95 drops from 4 to 2 between N=5 and N=7, then stays constant at 2 for N=9 and N=11. The number of active features falls from 9 to 8 at N greater than or equal to 7 because ph03 (the phase angle of the off-diagonal element) drops below the std=1e-6 threshold; under Markovian Z-dephasing the off-diagonal phase becomes effectively frozen.

PC1 has Purity as its dominant proxy at every tested N, but the PC1 loadings actually have similar magnitudes on Purity, von Neumann entropy, and several Bell fidelities. The reason PC1 covers so much variance at N=11 (74 percent) is that under post-ESD classical decoherence, Purity, SvN, and the Bell fidelities all evolve as deterministic functions of a single underlying mixture parameter, so they collapse onto one principal component.

PC2 at N=11 chain is loaded primarily on `psi_plus` (-0.69), `concurrence` (+0.51), and `phi_plus` (+0.38), and explains 22.7 percent of the variance. This is the "pre-ESD signature" axis: it captures the brief window between t=0 and the ESD time during which the Bell pair is still entangled and the Bell fidelities can move independently of the purity. This PC contributes meaningfully to coverage even though it lives in only a small fraction of the trajectory's wall-clock duration, because the variance it captures during the pre-ESD window is large.

Together, PC1 and PC2 cover 97.2 percent of the chain center_bell N=11 variance. PC3 adds another 1.8 percent. The cockpit's "3-observable" claim is satisfied with room to spare.

---

## 4. Star results

| N  | n_active | n95 | 3-PC coverage | PC1 variance | PC1 best proxy |
|----|----------|-----|---------------|--------------|----------------|
| 5  | 9        | 4   | 91.9%         | 51%          | Purity         |
| 7  | 9        | 4   | 93.0%         | 52%          | Purity         |
| 9  | 9        | 4   | 94.8%         | 54%          | Purity         |
| 11 | 9        | 3   | 95.3%         | 58%          | Purity         |

The star topology shows a flatter scaling profile than the chain. Coverage rises gently from 91.9 to 95.3 percent. n95 stays at 4 for N=5, 7, 9 and drops to 3 only at N=11. All 9 features remain active across every tested N: unlike the chain, the star center_bell trajectory keeps its phase angle ph03 above the std=1e-6 threshold. This is a topological consequence of the star's hub-and-spoke connectivity: the central hub qubit redistributes Hamiltonian-driven dynamics between Bell pair and spectator leaves, which keeps the off-diagonal phase oscillating long enough to register variance.

PC1 in star topology stays close to 50 to 58 percent of total variance across all N, indicating that the dominant axis is less dominant than in the chain case. The star trajectory is more genuinely multi-dimensional because the longer ESD time leaves more room for independent variation of the Bell fidelities and the concurrence.

---

## 5. The 8-feature reduction at chain N greater than or equal to 7

At chain N greater than or equal to 7 the analysis script drops one feature (ph03) before standardization because its variance falls below 1e-6. This is reported in the `n_active` column above as 8 instead of 9.

The empirical observation is that ph03 stays at exactly zero throughout the trajectory for chain center_bell pairs at N greater than or equal to 7 (numerically the standard deviation is around 6e-18, which is floating-point noise around exact zero). The mechanism is presumably a symmetry of the initial state preserved by the combined Hamiltonian and dephasing dynamics: Bell+ has a real off-diagonal element `rho[0,3] = 0.5`, and the combined Heisenberg-plus-Z-dephasing evolution appears to preserve this realness for chain initial states with `|+>^(N-2)` tensor factors. A formal proof of this symmetry is not attempted here; the empirical observation is sufficient for the gate to drop the feature correctly.

This is not a failure of the feature set, it is a correct identification by the gate that ph03 carries no usable PCA information for these specific configurations. In the star case, where the Hamiltonian dynamics are richer due to the hub coupling, ph03 stays above threshold and remains active.

A future iteration of the cockpit feature set could replace ph03 with something more informative for chain topologies, for example the magnitude of `rho[0,3]` (which is non-trivially decaying) instead of its phase. This is logged as a follow-up question, not a problem with the current result.

---

## 6. The far_edge control pair and what it actually showed

The far_edge pair was included in the experiment as an expected-trivial control class: two `|+>` qubits at the chain boundary, far from the central Bell pair, expected to show degenerate dynamics dominated by local dephasing without quantum-information content. The expectation was that the sanity gates would catch and drop these configurations as trivial.

What actually happened: at N greater than or equal to 7, the chain far_edge pair has concurrence exactly zero throughout the trajectory (max = 0, std = 2.5e-22, which is floating point noise around zero) and ph03 exactly zero (std = 6.6e-18, same situation). However, the other features (purity, von Neumann entropy, Bell fidelities, psi_norm) do show variation, because the two boundary `|+>` qubits undergo their own local dephasing dynamics. Purity ranges from 1.0 to 0.256 over the trajectory, exactly the same range as the center_bell pair, because the local decoherence of `|+>|+>` is structurally identical to the post-ESD phase of a Bell pair.

This means: 7 of 9 features are active, the purity range is large, Gates 2 and 3 both pass. Gate 1 does not apply to far_edge. The PCA on this trajectory gives n95=1 and 3-PC coverage above 99 percent, but the result is not informative about the cockpit framework's ability to monitor entangled observers; it is just confirming that classical decoherence is one-dimensional.

The far_edge pair therefore reports as "analyzed" rather than "dropped", but its high coverage number is not a contribution to the headline scaling result. The cockpit framework's relevant scope is the entangled observer class (center_bell), not arbitrary pairs in the system. Reporting the far_edge numbers in the same table as center_bell would be misleading; this section exists to make that explicit.

A future iteration of the sanity gates could extend Gate 1 to apply to all pairs, or add a new gate that requires concurrence variance to be non-trivial as a precondition for cockpit-relevance. This is logged as a pipeline improvement for the next experiment.

---

## 7. Adjacent pair (informational)

| N  | Chain adjacent 3-PC | Chain adjacent n95 | Chain PC1 proxy | Star center_leaf 3-PC | Star n95 | Star PC1 proxy |
|----|---------------------|--------------------|-----------------|-----------------------|----------|----------------|
| 5  | 89.6%               | 4                  | Purity          | 94.0%                 | 4        | Purity         |
| 7  | 95.1%               | 3                  | Purity          | 94.6%                 | 4        | Psi-norm       |
| 9  | 91.9%               | 4                  | Psi-norm        | 97.2%                 | 3        | Psi-norm       |
| 11 | 91.6%               | 4                  | Psi-norm        | 95.4%                 | 3        | Psi-norm       |

The adjacent pair is initially separable (one `|+>` qubit and one Bell-pair qubit, initial purity 0.5, initial concurrence 0) but is Hamiltonian-coupled to the Bell pair. Its trajectory therefore picks up dynamics indirectly. The 3-PC coverage stays in the 89 to 97 percent range across all N for both topologies, comfortably above the 85 percent threshold.

The interesting feature here is the PC1 proxy transition from Purity to Psi-norm at N=9 for chain (and at N=7 for star). At small N, the dominant axis of variation is bulk mixture (Purity). At larger N, it shifts to coherence magnitude (Psi-norm), reflecting that as the chain or star grows the adjacent pair's dynamics become coherence-driven rather than mixture-driven. This is a qualitative transition, not just a quantitative shift, and it would be worth following up in a separate experiment that asks whether the adjacent pair's PC1 identity transition has a precise N at which it occurs (and whether it corresponds to a topological or spectral threshold).

---

## 8. Cavity reframing: the same result in optical language

Everything in this document so far has been written in the quantum-information vocabulary that the original COCKPIT_UNIVERSALITY paper used: concurrence, entanglement, Bell pair lifetime, monogamy of entanglement. That language is correct, but it is not the language the rest of this repository uses. Since early April 2026 the canonical framing has been the optical cavity reframing introduced in [OPTICAL_CAVITY_ANALYSIS](OPTICAL_CAVITY_ANALYSIS.md) and developed in [VEFFECT_CAVITY_MODES](VEFFECT_CAVITY_MODES.md), [CAVITY_MODE_LOCALIZATION](CAVITY_MODE_LOCALIZATION.md), and the [Absorption Theorem](ABSORPTION_THEOREM_DISCOVERY.md). In that language, the cockpit scaling result is not a new finding so much as a direct empirical application of theorems that already exist in the repo.

This section translates the result into cavity language and shows the connections.

### Translation table

| Quantum-information language | Optical cavity language |
|------------------------------|-------------------------|
| Bell pair initial state on `(c1, c2)` | Coherent input mode pumped into the cavity at qubits `c1, c2` |
| Concurrence of the Bell pair | Magnitude of the cavity's coherent (light) component |
| Entanglement Sudden Death | Absorption of the coherent mode by the illumination `gamma` |
| ESD time | Mode absorption time, i.e. cavity finesse for this specific mode |
| Z-dephasing rate gamma | External illumination intensity per site |
| Heisenberg coupling J | Internal cavity coupling between elements |
| Markovian limit (N greater than or equal to 7) | Bath large enough that absorbed light does not return |
| Entanglement revival (chain N=5) | Non-Markovian re-emission from a finite cavity |
| Monogamy of entanglement | Aperture distribution of light over multiple collection elements |
| n95 (effective dimensionality) | Effective number of cavity modes the trajectory traverses |
| Purity as PC1 proxy | Mode purity as the dominant trajectory axis |
| post-ESD classical phase | Post-absorption classical diffusion regime |
| pre-ESD quantum phase | Pre-absorption coherent mode lifetime |

### The Bell pair as cavity input

The Bell+ state on the central pair has a clean Pauli decomposition: `rho_Bell+ = (1/4)(II + ZZ + XX - YY)`. Three of the four Pauli strings are weight-zero or pure-Z (`II`, `ZZ`), one weights as `n_XY = 2` per string (`XX` and `YY` both contain two transverse Pauli factors). In the cavity language, the Bell pair is half "structure" (the `II + ZZ` content, immune to illumination) and half "light" (the `XX - YY` content, fully exposed to absorption). The entangled part of the Bell pair is precisely the light-bearing part: concurrence reads exactly the same Pauli content that the Absorption Theorem says will be absorbed.

This is why ESD looks like sudden absorption rather than smooth decay. The light-bearing components `XX` and `YY` are not eigenmodes of the Lindbladian; they get mixed by the Hamiltonian into multiple cavity modes, each of which has its own `n_XY` and absorbs at rate `2*gamma*n_XY`. The Bell pair concurrence is a non-linear function of how much of the original light content survives across all those modes simultaneously. When enough of the light is absorbed, the concurrence drops to zero in finite time, even though the underlying mode magnitudes are still decaying smoothly.

### The chain vs star asymmetry as an aperture effect

The Absorption Theorem says: every Liouvillian eigenmode has decay rate `Re(lambda) = -2*gamma*<n_XY>`. The implication for the Bell pair is that the absorption time of its light-bearing components depends on how the Hamiltonian distributes the `XX, YY` content across the eigenmodes. If the relevant modes have high `<n_XY>` (concentrated light), absorption is fast. If the modes have low `<n_XY>` (distributed light), absorption is slow.

**Chain cavity.** A linear arrangement, `N - 1` bonds. The Hamiltonian mixes the Bell pair light into modes that span the chain, but the light density per mode is set by the local geometry around the Bell pair position and is approximately N-independent for a center pair (the boundary qubits are too far away to matter). The dominant Bell-pair-bearing modes therefore have approximately constant `<n_XY>` regardless of N. Absorption rate is constant. Bell pair lifetime (ESD time) is constant at t ~ 1.

**Star cavity.** A hub with `N - 1` leaves, all bonds going through the same central qubit. The Bell pair sits on two leaves, but because every leaf is coupled to the hub, the Hamiltonian distributes the Bell pair light across **all leaves** through the central node. With more leaves (larger N), the same total light is spread across more cavity modes. Each individual mode has lower `<n_XY>`, lower absorption rate, longer lifetime. The Bell pair lifetime grows with N because the cavity has a larger collection aperture for distributing the light over more modes.

This is what the [V-Effect cavity reframing](VEFFECT_CAVITY_MODES.md) already showed at N=5 with a different observable: chain has 112 distinct frequencies but Q_max = 72.4, while star has only 42 frequencies but Q_max = 100.0. Star has fewer modes but each one lives longer. Today's cockpit scaling result is the time-domain version of that frequency-domain finding: chain has fast absorption per mode and short ESD, star has slow absorption per mode and long ESD.

### What `n95` actually counts

In the quantum-information language, `n95` is the number of principal components needed to capture 95 percent of the trajectory variance. In the cavity language, it is the number of distinct cavity modes whose evolution the trajectory traces out before absorption renders them indistinguishable.

A low `n95` (chain post-ESD, n95 = 2) means: after absorption, the cavity is in a single dominant classical diffusion mode (PC1, dominated by purity decay) plus a small residual coherent component captured by PC2. The trajectory has effectively collapsed to two dimensions because all the originally distinct light-bearing modes have been absorbed and the remaining dynamics is one-dimensional classical mixing.

A higher `n95` (star, n95 = 3 to 4) means: the absorption is slower, so multiple light-bearing modes remain distinguishable for a substantial fraction of the trajectory. PCA picks them up as independent principal components.

The cockpit framework's "first 3 PCs cover 88 to 96 percent" claim from COCKPIT_UNIVERSALITY is, in this language, the statement that any cavity dominated by absorption settles to a low-mode regime within a few absorption times. The 3-observable cockpit is sufficient because the cavity dynamics is structurally sparse, not because the observables are particularly clever.

### Purity as the natural cavity coordinate

Purity = `Tr(rho^2)` measures how much of the cavity mass is concentrated in a single quantum state versus distributed across many. In cavity language, this is **mode purity**: the inverse participation ratio of the cavity's instantaneous state across the basis of eigenmodes. Purity equals 1 when the cavity is in a single coherent mode. It drops as the mode content spreads out under absorption and Hamiltonian mixing.

The reason Purity is the universal PC1 proxy across topologies and across N is that it directly reads the cavity's mode purity, which is the natural coordinate for any cavity description. The cockpit framework is, at its core, a Purity-driven monitor, with the other observables (Bell fidelities, von Neumann entropy, concurrence) serving as orthogonal corrections that pick up the residual dynamics not captured by the mode purity axis.

### Why this matters for what comes next

Three immediate consequences of taking the cavity reframing seriously:

1. **The Absorption Theorem is the natural predictor of ESD time, not an external reference.** Instead of measuring ESD empirically, we should be able to compute it from the spectrum of the Liouvillian directly: identify the slowest-absorbing mode that has support on the Bell pair light components, and read its absorption rate as `2*gamma*<n_XY>`. This is the analytical follow-up the experiment has been waiting for, and it lives entirely within the existing repo theorems.

2. **The topology lever is an aperture lever.** The factor-of-4 lifetime improvement of star over chain at N=11 is not "monogamy magic". It is the same kind of improvement an optical system gets by widening its collection aperture. Topologies that distribute light over many modes win. This unifies the cockpit scaling result with the existing sacrifice zone work (which is also about distributing absorption optimally, just along a different axis: per-site instead of per-mode).

3. **The Bell pair is not the only useful input.** The cockpit framework has been using Bell+ as the canonical entangled observer, but the cavity reframing makes it natural to ask which input states maximize the cavity's response. A pure-light input (something with `<n_XY>` close to N) would absorb fastest. A pure-structure input (something close to `II + ZZ`) would never absorb at all. The cockpit's effective dimensionality is a function of the input choice, not just of the system. This connects directly to [WHAT_QUBITS_EXPERIENCE](../hypotheses/WHAT_QUBITS_EXPERIENCE.md) and [GAMMA_IS_LIGHT](../hypotheses/GAMMA_IS_LIGHT.md).

The honest summary of this section: the cockpit scaling result of N=5 to N=11 is, in retrospect, a numerical confirmation of two theorems already proved in the repo (the Absorption Theorem from April 4, and the V-Effect cavity geometry from April 4). The data was new today; the explanation was waiting in the repo since three days ago. The next iteration of cockpit work should start from the cavity language and use the Absorption Theorem as the predictive tool, with quantum-information observables (concurrence, n95) as derived measurements rather than primary objects.

---

## 9. Verdict

**The cockpit framework scales for the entangled observer class, mediated by Entanglement Sudden Death.**

For Bell pair observers in Heisenberg systems under uniform local Z-dephasing, the 3-observable cockpit framework continues to capture more than 90 percent of the trajectory variance across N=5 to N=11 in both chain and star topologies. The effective dimensionality n95 decreases for chain (4 to 2) and decreases more gently for star (4 to 3), but in both cases stays at or below the original small-N values rather than growing linearly with N as the COCKPIT_UNIVERSALITY extrapolation tentatively suggested.

The decrease in n95 is not a "magical simplification" of the dynamics. It is a consequence of Entanglement Sudden Death: once the Bell pair's concurrence collapses to zero (which happens at t ~ 1 for chain regardless of N, and at t ~ 0.5 to 3.9 for star depending on N), the rest of the trajectory is classical decoherence on a separable mixed state, which is structurally low-dimensional. The cockpit's "first 3 PCs" are sufficient because the trajectory itself only has 2 to 4 effective dimensions. The framework's claim ("3 observables suffice") was conservative for the post-ESD phase; the new finding is that this remains true even when the pre-ESD phase contributes a meaningful but small variance fraction.

Purity is the dominant PC1 proxy in every analyzed configuration. This is consistent with the original COCKPIT_UNIVERSALITY result (which found Purity dominant for chain at N greater than or equal to 4) and extends it to N=11 and to star topology.

**What is confirmed:**
- 3-observable coverage stays above 90 percent for the entangled observer class up to N=11
- Purity is the dominant cockpit observable across topologies and system sizes
- The framework is consistent with known physics (ESD, monogamy of entanglement) rather than dependent on accidental low-dimensional embeddings

**What is qualified:**
- The framework's effective dimensionality at large N is set by ESD timing, which is topology-dependent
- The cockpit's "richness" comes increasingly from the post-ESD classical phase as N grows; the pre-ESD quantum phase contributes to PC2 but its time-share decreases
- The n95 numbers reported here are conservative for the entangled observer class only; arbitrary pairs in the system can show artificially low n95 due to one-dimensional classical decoherence (see Section 6)

---

## 10. Limitations

1. **Heisenberg interactions and Z-dephasing only.** All results assume the standard Heisenberg coupling (XX+YY+ZZ) and uniform local Z-dephasing. Other coupling schemes (XX-only, anisotropic, long-range) and other noise models (depolarizing, amplitude damping, non-Markovian) are not tested. The COCKPIT_UNIVERSALITY baseline included depolarizing noise at N=2-4 and showed similar dimensionality, but extending depolarizing tests to N=11 is a separate task.

2. **N=11 is the largest tested.** The C# matrix-free propagator can handle N=15, but the runs were not executed for this experiment because the trend stabilized clearly at N=11. A future experiment could verify that the chain n95=2 plateau continues to N=15, but the prediction is that it will, since the ESD-driven mechanism does not depend on system size beyond the threshold where the bath becomes effectively Markovian.

3. **Three pair types per configuration.** Only the center_bell, adjacent, and far_edge (or center_leaf, far_leaf for star) pairs are extracted. The full pair-distance scan that is available at N=5 in COCKPIT_UNIVERSALITY (all 10 pairs) is not reproduced here. This is sufficient to answer the scaling question, but a comprehensive distance-resolved scaling map would require extracting more pairs per N.

4. **The far_edge control class is not informative for the cockpit claim.** As discussed in Section 6, the sanity gates allow far_edge pairs to be analyzed but their reported coverage is not a measurement of the cockpit framework's relevant scope. Future iterations of the gates should make this distinction explicit at the pipeline level rather than at the documentation level.

5. **The adjacent pair PC1 proxy transition (Section 7) is observed but not characterized.** It would be worth a separate experiment to find the exact N at which the transition occurs and whether it corresponds to a topological or spectral feature of the underlying Heisenberg system.

6. **The ph03 freezing under chain Z-dephasing (Section 5) is observed but not formally proved.** The phase angle of the off-diagonal element stays at exact zero throughout the trajectory for chain center_bell pairs at N greater than or equal to 7, presumably due to a symmetry of the initial state preserved by the combined dynamics. A formal derivation is left as an open question.

---

## 11. Cross-validation against direct spectral evolution (April 7, 2026)

The N=5 chain center_bell trajectory was independently reproduced by direct spectral evolution of the full Liouvillian, as a sanity check on the C# matrix-free propagator and as a concrete application of the Absorption Theorem (Section 8). The script [path_d_bell_pair_absorption.py](../simulations/path_d_bell_pair_absorption.py) builds the Heisenberg Hamiltonian, constructs the 1024 by 1024 Liouvillian (column-stacking vectorisation, so each Pauli string becomes one basis vector in the operator space), and diagonalises it. Every right eigenoperator is decomposed in the Pauli basis to obtain the weighted average `<n_XY>` per mode. The initial state `|Phi+><Phi+|_{2,3}` tensor `|+><+|^{otimes 3}` is then projected onto the eigenoperators and evolved exactly to `rho(t)` at 301 time points between t=0 and t=3. Reduced 2-qubit density matrices and Wootters concurrence are computed from the result.

**Theorem verification.** The Absorption Theorem `Re(lambda) = -2 * gamma * <n_XY>` holds for all 1024 Liouvillian modes of the N=5 chain to a maximum deviation of 2.4e-14, consistent with floating-point precision. This is a redundant check (the theorem has already been verified across 1343 modes in [ABSORPTION_THEOREM_DISCOVERY](ABSORPTION_THEOREM_DISCOVERY.md)), but it confirms that the spectral pipeline used here is consistent with the existing proof.

**Trajectory match.** The spectral concurrence trajectory matches the empirical CSV from `cockpit_scaling_v2_N5_chain.csv` to a maximum deviation of 0.001144 (RMS 0.000353) across 31 sampled time points, well within CSV write precision. The match includes the smooth Markovian decay from t=0 to t=0.9, the entanglement sudden death at t=0.9, the entanglement revival to a peak of 0.158 at t=2.1, and the final death at t=2.4. Purity matches to the same precision.

![Path D cross-validation: empirical cockpit_scaling_v2 trajectory (red circles) vs direct spectral evolution (blue curve) for the Bell pair on qubits 2,3 of the N=5 Heisenberg chain at gamma=0.05. Concurrence (top) and purity (bottom) are reproduced to within 0.001 across all sampled time points, including the entanglement sudden death at t=0.9 and the revival peak of 0.158 at t=2.1.](../simulations/results/path_d_comparison.png)

**Open analytical question.** The spectrum reveals a clean structural picture for the Bell pair as input. There are exactly 50 Liouvillian modes sitting at `Re(lambda) = -0.2`, the n_XY=2 absorption plateau predicted by the theorem. All 50 have `<n_XY>` exactly equal to 2.000000 to floating-point precision, confirming they form a sharply defined sector rather than a smeared band. These 50 modes together carry 75 percent of the Bell pair's projection onto the broader n_XY=2 sector (defined as modes with `<n_XY>` between 1.5 and 2.5), and within that 50-mode plateau the distribution is heavily skewed: just 6 modes carry 90 percent of the plateau-internal weight, with the top two alone accounting for about 70 percent. Despite this clean spectral picture, predicting the ESD time t=0.9 from the spectrum alone is not straightforward: taking the slowest plateau mode and computing `-log(threshold) / rate` gives t_ESD around 88, off by a factor of about 100. The actual ESD time is governed by destructive interference among the dominant plateau modes through the Wootters concurrence functional, which is intrinsically nonlinear in the density matrix. A clean analytical formula for ESD time as a function of the spectrum (in the spirit of the Absorption Theorem itself) does not yet exist for entangled observers and remains an open question.

**Robustness note.** numpy.linalg.eig chooses an arbitrary basis within each degenerate eigenvalue subspace, depending on the BLAS backend, so the per-mode overlap distribution among the 50 plateau modes is not unique across machines. The concurrence trajectory is identical regardless, because Wootters concurrence is invariant under unitary rotations within degenerate subspaces. The cockpit framework is robust because it does not depend on identifying any individual mode as special.

**What this cross-validation establishes:**
1. The cockpit_scaling_v2 trajectories are numerically equivalent to the textbook Lindblad spectral solution. The C# matrix-free propagator and the spectral pipeline produce the same dynamics with no hidden numerical drift in either.
2. The Absorption Theorem alone correctly identifies the absorption rate of every Liouvillian mode, but it is not yet sufficient to predict the ESD time of an entangled observer. ESD is a threshold effect on a nonlinear functional of the density matrix, and that nonlinearity is not yet folded into a closed-form prediction.

---

## 12. References

- [COCKPIT_UNIVERSALITY](COCKPIT_UNIVERSALITY.md) -- the N=2-5 baseline result, the framework definition
- [cockpit_scaling_analysis_v1.py](../simulations/cockpit_scaling_analysis_v1.py) -- the Python analysis script with sanity gates and dropping logic
- [Program.cs cockpit mode](../compute/RCPsiSquared.Propagate/Program.cs) -- the C# cockpit dispatch and trajectory generator
- [DensityMatrixTools.cs](../compute/RCPsiSquared.Propagate/DensityMatrixTools.cs) -- BellFidelity, Ph03, ExtractCockpitFeatures helpers
- [cockpit_scaling_v2_results.txt](../simulations/results/cockpit_scaling_v2/cockpit_scaling_v2_results.txt) -- the full numerical output
- [cockpit_scaling_v2_results.json](../simulations/results/cockpit_scaling_v2/cockpit_scaling_v2_results.json) -- per-configuration JSON dump for re-analysis
- [cockpit_scaling_v2_curve.png](../simulations/results/cockpit_scaling_v2/cockpit_scaling_v2_curve.png) -- the four-panel scaling figure embedded above
- [path_d_bell_pair_absorption.py](../simulations/path_d_bell_pair_absorption.py) -- the spectral cross-validation script described in Section 11
- [path_d_comparison.png](../simulations/results/path_d_comparison.png) -- visual overlay of empirical vs spectral concurrence and purity trajectories
- [PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) -- the Lindblad rate structure that directly determines ESD timing in this experiment via `Re(lambda) = -2*gamma*<n_XY>`
- [ABSORPTION_THEOREM_DISCOVERY](ABSORPTION_THEOREM_DISCOVERY.md) -- the empirical discovery of the absorption identity, the source of the `<n_XY>` interpretation used in Section 8
- [VEFFECT_CAVITY_MODES](VEFFECT_CAVITY_MODES.md) -- the topology Q-factor table (chain N=5: Q_max=72.4, star N=5: Q_max=100.0) that this experiment confirms in the time domain
- [OPTICAL_CAVITY_ANALYSIS](OPTICAL_CAVITY_ANALYSIS.md) -- the original Fabry-Perot reframing that introduced the cavity language used in Section 8
- [GAMMA_IS_LIGHT](../hypotheses/GAMMA_IS_LIGHT.md) -- the hypothesis that gamma is external illumination, central to the Section 8 reading

---
