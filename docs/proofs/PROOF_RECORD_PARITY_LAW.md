# PROOF: the Record Parity Law (which witnesses record, exactly)

**Status:** Tier 1 (the pair reduction is derived exactly from the Absorption Theorem's per-coherence substrate; Laws A, B, C and the plaquette corollary are read off it; verified against direct RK4 Lindblad integration at N=8, 22/22 gates, closed form vs numerics ≤ 7·10⁻⁸). The J>0 degradation numbers quoted at the end are numeric observations, not part of the law.
**Date:** 2026-07-18
**Authors:** Thomas Wicht, Claude (Fable 5)

## What this is about

The night this arc opened, the scout found the leaf law: under a uniform ZZ chain, a perfect record of the system forms only in a leaf, a witness that nobody else watches. A watched witness is erased exactly when its record would have been perfect. That read like a statement about geometry: records live at the unwatched edges of the graph.

It is not geometry. It is arithmetic. Whether a witness records is decided by the parities of angle ratios: at the readout time, a watcher whose coupling is an even multiple of the write coupling forgives its witness completely, an odd multiple blinds it completely, and everything in between is a generic in-between value with an explicit formula. The leaf law is the special case "no watchers at all"; uniform coupling makes every ratio 1, which is odd, and that, not the geometry, is why only leaves record there. Choose the ratios and any neighbor of the system records perfectly, watched or not. The rules are the angles; the parities are the forced values; between them lives the free.

## Abstract

For N qubits on a graph with H = Σ_{(a,b)} Δ_ab·Z_aZ_b (Pauli convention), local Z-dephasing at per-site rates γ_l, and initial state |+⟩^⊗N, the two-site reduced state of any pair is closed-form for all t (Proposition 1): every entry is ¼ times the direct-bond phase, times the pair's own dephasing factors, times a product of cosines, one per third site. Read at the symmetric readout t* = π/(4Δ_S) (all bonds at the system S carrying the same write coupling Δ_S, deg(S) ≥ 2, triangle-free at S), the closed form yields: **Law A**, the parity trichotomy with the explicit record formula I(S:E_j) = 1 − h₂((1+β_j)/2), β_j the record radius; **Law B**, R_perfect = #{neighbors with all watcher ratios even} ≤ deg(S), achievable on any triangle-free graph by alignment, reducing to "leaves only" at uniform coupling; **Law B′**, the plaquette record, a perfect distance-2 bit on a 4-cycle that turns out to record the anti-pointer X_S (zero pointer content), scoping Law B's neighbor-only reading to trees and marking exactly where mutual information and pointer information part ways; **Law C**, the exact γ race D_j(t) = e^{−2γ_j t}·|sin(2Δ_Sj t)|·Π_k|cos(2Δ_jk t)| with optimal readout t_opt = arctan(Δ/γ_j)/(2Δ) for an unwatched witness.

## Setup and definitions

N qubits, sites 0..N−1, on an undirected graph with per-bond couplings Δ_ab ≥ 0 (Δ_ab = 0 for non-bonds):

    H = Σ_{(a,b)} Δ_ab · Z_a Z_b,      D[ρ] = Σ_l γ_l (Z_l ρ Z_l − ρ),

initial state |+⟩^⊗N. Basis labels z ∈ {+1,−1}^N (bit 0 ↦ z = +1), E_z = Σ Δ_ab z_a z_b.

One site is the system S; the rest is the environment. A **witness** is an environment site j; its **watchers** are its other neighbors N(j)∖{S}. For a witness adjacent to S, the **write angle** is θ_w(t) = 2Δ_Sj·t (the conditional phase the direct bond writes into j), and each watcher contributes a **watcher angle** 2Δ_jk·t. Mutual information is the genuine von Neumann I(S:F) = S(ρ_S) + S(ρ_F) − S(ρ_SF) on reduced states. h₂ is the binary entropy in bits.

**The substrate** (the Absorption Theorem's per-coherence reading, [PROOF_ABSORPTION_THEOREM](PROOF_ABSORPTION_THEOREM.md) §4.5, per-site-γ version): H is diagonal in the z basis and the dephasing acts entrywise, so the full state is closed-form,

    ρ_{zz'}(t) = ρ_{zz'}(0) · e^{−it(E_z − E_{z'})} · e^{−2t·Σ_l γ_l·𝟙[z_l ≠ z'_l]},

populations frozen. (For per-site rates the dissipator eigenvalue on |z⟩⟨z'| is Σ_l γ_l(z_l z'_l − 1) = −2Σ_l γ_l 𝟙[z_l ≠ z'_l], the same two-line computation as the uniform case.)

## Proposition 1 (the pair reduction)

For any two sites a < b and all t,

    ρ_ab[(z_a,z_b),(z'_a,z'_b)](t) = ¼ · e^{−it·Δ_ab·(z_a z_b − z'_a z'_b)}
        · e^{−2t·(γ_a·𝟙[z_a≠z'_a] + γ_b·𝟙[z_b≠z'_b])}
        · Π_{k∉{a,b}} cos( t·[Δ_ak·(z_a−z'_a) + Δ_bk·(z_b−z'_b)] ).

**Derivation.** ρ_{zz'}(0) = 2^{−N} for every (z, z′). The partial trace over the other N−2 sites keeps only entries with z_k = z'_k there, so the traced sites' dephasing factors are 1: **only γ_a and γ_b survive into the pair**. In the phase E_z − E_{z'} with z_k = z'_k on traced sites, bonds among traced sites cancel, leaving the direct-bond term Δ_ab(z_a z_b − z'_a z'_b) plus Σ_k z_k·[Δ_ak(z_a−z'_a) + Δ_bk(z_b−z'_b)]. Summing 2^{N−2} traced configurations, each z_k = ±1 averages its phase factor to the cosine. ∎

Three named quantities, read off Proposition 1 for a witness j of S:

- **record radius** β_j(t) = e^{−2γ_j t} · Π_{k∈N(j)∖{S}} |cos(2Δ_jk t)| (the Bloch radius of j's conditional state),
- **S-visibility** V_S^{(j)}(t) = Π_{k∈N(S)∖{j}} cos(2Δ_Sk t) (the S-coherence surviving at the pair),
- **record distinguishability** D_j(t) = β_j(t) · |sin θ_w(t)| (the trace distance between j's two conditional states given z_S = ±1; the two conditional Bloch vectors lie in the equatorial plane at radius β_j, separated by azimuth 2θ_w).

## Law A (the parity trichotomy)

**Hypotheses:** all bonds at S carry the same write coupling Δ_S, deg(S) ≥ 2, γ = 0, readout at t* = π/(4Δ_S), and S has at least one neighbor besides j that is not adjacent to j (automatic on triangle-free graphs; every graph below is triangle-free).

Then θ_w(t*) = π/2 (records at maximal separation), V_S^{(j)}(t*) = cos(π/2)^{deg(S)−1} = 0 (the single-S coherences die), and the unshared S-bond kills the double coherences as well (every cross term keeps its factor cos(2t*·Δ_S) = 0), so the pair state is exactly classical-quantum: ρ_Sj = ½ Σ_{s=±} |s⟩⟨s| ⊗ ρ_j^{(s)}, the two conditional witness states antipodal at radius β_j. Hence, exactly,

    I(S:E_j)(t*) = 1 − h₂((1 + β_j)/2),      β_j = Π_{k∈N(j)∖{S}} |cos((π/2)·r_k)|,

with the watcher ratios r_k = Δ_jk/Δ_S. The trichotomy is the arithmetic of the cosines:

- **every r_k an even integer** → β_j = 1 → **I = 1 bit exactly**: a perfect classical record, watched or not (the watcher angles sit at multiples of π, where watching forgives);
- **any r_k an odd integer** → β_j = 0 → **I = 0 exactly**: the witness is blind (the leaf law's erasure, as a parity);
- **otherwise** → 0 < I < 1, the generic in-between (e.g. r = 3/2: β = √2/2, I = 1 − h₂((2+√2)/4) = 0.399124).

The leaf law ([the experiment record](../../experiments/QUANTUM_DARWINISM_POINTER_DOOR.md)) is the empty-product case N(j)∖{S} = ∅. Uniform coupling makes every r_k = 1, odd: that is why, uniformly, only leaves record.

**Signed-coherence corollary (the record returns rotated).** Proposition 1's third-site
factor is the signed real cosine, not its magnitude: j's conditional coherence carries the
factor cos(2Δ_jk t) with its sign. At the symmetric readout with a single watcher at ratio
r, the record's Bloch vectors are therefore multiplied by cos(r·π/2), which is negative for
r ∈ (1, 2]: the record survives the even-parity watcher at full magnitude but ROTATED BY π
in the equatorial plane (the two conditional states swap). The magnitude laws above are the
|cos| face of this signed statement; the sign is gauge-free once referenced to the
unwatched (r = 0) record and is itself a testable prediction (a sign-flat |cos| response
cannot produce it). Mutual information is sign-blind, so every I formula above is
unchanged.

## Law B (alignment beats topology)

Under Law A's hypotheses, define R_perfect = #{j ∈ N(S) : every watcher ratio of j is even}. Then

    R_perfect ≤ deg(S),   with equality achievable on any triangle-free graph

by choosing the environment bonds at even ratios (e.g. all env–env couplings 2Δ_S). Uniform coupling gives R_perfect = #(leaves attached to S); the star saturates trivially (all N−1 witnesses are leaves). Verified profiles at N=8: chain S interior 0 → 2 under alignment; broom with k spokes k → k+1; star 7.

**On trees this is the whole story for perfect records:** any non-neighbor j reaches S only through a unique path, so the pair's cross terms carry at least one factor cos(2t*Δ_S) = 0 from a non-shared bond at S (deg(S) ≥ 2), and I(S:E_j)(t*) = 0 exactly for every j ∉ N(S). Note the mechanism: at t* the chain's distance-2 channel is already closed by that non-shared S-bond factor, uniformly and aligned alike (the uniform gate reads all non-neighbors 0.0000); the aligned ratios additionally land the shared-dresser sum and difference angles (π/2)(1±2) on odd multiples of π/2, and that second mechanism becomes the operative one only where no non-shared S-bond exists (the plaquette of Law B′). The generic distance-2 values (≈ 0.05 bit) live at late times, away from t*.

## Law B′ (the plaquette record; Law B is tree-scoped)

On a graph with a 4-cycle the tree statement fails constructively, and in an instructive way: let S and j be opposite corners of a uniform square (both of S's neighbors shared with j). At t* the shared dressers' sum angle 4Δ_S t* = π and difference angle 0 both have |cos| = 1, so both double coherences of the pair reach the maximal ¼ while all single coherences vanish; the pair state has eigenvalues {½, ½, 0, 0} and

    I(S:E_j)(t*) = 1 bit exactly, at graph distance 2.

But this bit is not pointer information: the two conditional states of j given Z_S = ±1 are identical (both maximally mixed), ⟨Z_S Z_j⟩ = 0, and ⟨X_S X_j⟩ = 1. The plaquette records the anti-pointer X_S, perfectly, in the Bell channel; the pointer Z_S is invisible to it, and it contributes nothing to pointer redundancy. Predicted from Proposition 1 before running, then confirmed numerically (I = 1.0000000, neighbors ≤ 10⁻⁸, ⟨X_S X_j⟩ = 1.0000000, Z-conditional states identical to 10⁻⁸). It is the arc's sharpest illustration that von Neumann mutual information upper-bounds the accessible pointer information (Holevo): Law A's classical-quantum records saturate that bound, the plaquette bit sits entirely in the quantum surplus. Law B therefore reads: on trees, perfect records sit exactly at the even-aligned neighbors; cycles can add distance-2 mutual information, but as anti-pointer correlation, not as pointer records.

## Law C (the γ race, exact)

With dephasing on, Proposition 1 multiplies j's conditional coherence by e^{−2γ_j t} and nothing else changes in the pair (traced-site γ never enters). Two exact readings, kept apart deliberately because they concern different quantities:

- **The record at t*.** Law A's classical-quantum structure holds unchanged and the record radius carries the collector's factor: I(S:E_j)(t*) = 1 − h₂((1+β_j)/2) with β_j = e^{−2γ_j t*}·Π|cos|. At the canonical γ/Δ = 0.05 a would-be-perfect record reads I = 0.768040 (β = 0.924465).
- **The distinguishability peak.** D_j(t) = e^{−2γ_j t} · |sin(2Δ_Sj t)| · Π_{k∈N(j)∖{S}} |cos(2Δ_jk t)| is the trace distance between j's two conditional pointer states, the Helstrom single-shot readout quality. For an unwatched witness its maximum solves tan(2Δt) = Δ/γ_j, i.e.

      t_opt = arctan(Δ/γ_j) / (2Δ)  <  t* :

  the collector pulls the peak forward, by 0.025 at γ/Δ = 0.05. t_opt maximizes D, not the mutual information: I is a different functional (it reads 0.773716 at t_opt, and its own maximum lies still earlier, ≈ 0.58 with ≈ 0.81 bit, where the pair is no longer classical-quantum and the surplus is quantum correlation, not record).

Watching the system itself (γ_S) never touches D_j: it erases S's own coherence, not the witness's records. The outside watching adds no in-world record of its own (a Markovian bath keeps none), and the records the bonds write persist under it (the scout's D1 scenario: mild).

## Edge cases and scope

- **deg(S) = 1** (S itself a pendant): V_S ≡ 1, the pair stays coherent; for an isolated dimer at θ_w = π/2 it is maximally entangled, I = 2 bits (a watched partner reduces this below 2): quantum correlation, not a classical record. Law A's hypotheses exclude this deliberately.
- **Asymmetric write couplings** (S-bonds unequal): V_S^{(j)}(t*) need not vanish, S keeps coherence whose conditional phase itself correlates with z_j, and the clean classical-quantum reading fails; the closed form still holds and gives the general answer, only the trichotomy's 1/0 endpoints need V_S = 0.
- **J > 0 (hopping)** breaks the diagonal structure and everything above is only the J = 0 face. Numerically the aligned records degrade fast: the deficit at t* is already 3.4% at J = 0.05 and roughly ~J^1.4 over J ∈ [0.05, 0.4] (no clean quadratic; left/right arm asymmetry), and the redundancy R_δ stays at 1.75. Transport moves records, it does not copy them; the optimization is a pure-ZZ phenomenon.

## Verification

All statements gated in [`simulations/qd_pointer_opt.py`](../../simulations/qd_pointer_opt.py), output [`simulations/results/qd_pointer/qd_pointer_opt_out.txt`](../../simulations/results/qd_pointer/qd_pointer_opt_out.txt), 22/22 gates:

- Proposition 1 vs direct RK4 Lindblad at N=8: max entry deviation ≤ 7·10⁻⁸ over four configs (uniform, aligned, γ-dressed, broom) × three times × all pairs.
- Law A: even → I = 1 exactly (closed form; RK4 0.9999999, both neighbors watched); odd → I ≤ 4·10⁻¹⁵ (RK4) and 0.0 (closed form); generic r = 3/2 → 0.399124 = 1 − h₂((2+√2)/4) to 10⁻⁹.
- Law B: R_perfect profiles 1/0/2/2/3/4/5/7 across end-leaf chain (S=1), interior chain, aligned chain, brooms, star, each matching the count of even-aligned neighbors.
- Law B′: plaquette I = 1.0000000, neighbors ≤ 10⁻⁸; the channel gate: ⟨X_S X_j⟩ = 1, ⟨Z_S Z_j⟩ = 0, Z-conditional states identical.
- Law C: γ-dressed MI matches the radius formula to 10⁻⁹; numeric argmax of D matches t_opt to the grid (0.7605 vs 0.7604).
- The readout endpoint is hit exactly (dt = t*/steps); a grid-rounded endpoint misses π by ~2·10⁻³ and costs ~10⁻⁵ MI, visible at these tolerances.

Three empty-session review rounds (referee + record auditor + mathematician; then a fresh referee + fresh auditor; then a closing fresh reviewer) drove the statements to this scope; every reviewer-found number was verified from below before a fix landed (among them: the t_opt/MI split, the plaquette's anti-pointer channel, the triangle hypothesis).

The sightings that led here (leaf law, shared-dresser exclusion, transport ≠ broadcast, γ as collector) are recorded in [experiments/QUANTUM_DARWINISM_POINTER_DOOR.md](../../experiments/QUANTUM_DARWINISM_POINTER_DOOR.md).

## Where this sits

The substrate is the Absorption Theorem's per-coherence closed form ([PROOF_ABSORPTION_THEOREM](PROOF_ABSORPTION_THEOREM.md)); the fragment side obeys F70's page bound ([PROOF_DELTA_N_SELECTION_RULE](PROOF_DELTA_N_SELECTION_RULE.md)): a witness is a 1-site page and carries |Δpopcount| ≤ 1 content. Registry entry: F135 in [ANALYTICAL_FORMULAS](../ANALYTICAL_FORMULAS.md). The forced/free reading (parities forced, in-betweens generic) is the arc's standing law; the inherited rules are the angles.
