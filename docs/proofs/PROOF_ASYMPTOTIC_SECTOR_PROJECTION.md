# Proof: Asymptotic Sector Projection Theorem

**Status:** Fully proven (all three steps analytic). Steps 2–3 additionally verified numerically.
**Date:** April 12, 2026
**Authors:** Thomas Wicht, Claude (chat)
**Context:** Identified during Track B of the three-values investigation (chat session 2026-04-12). Recorded here before the companion Track B experiment ([Sector Projection Formula](../../experiments/SECTOR_PROJECTION_FORMULA.md)) wrote it, so the experiment could cite the proof rather than re-derive it.

---

## Abstract

Under Heisenberg dynamics with local Z-dephasing at every site of any connected graph, the long-time state is decided entirely by where the initial state sits among the excitation sectors, not by any coherence between them. The excitation-sector populations p_w = Tr(P_w ρ₀) are exact constants of motion: the XX+YY hopping moves excitations within a sector without changing the total count, and the dephasing dissipator commutes with the sector projector. Every inter-sector coherence decays, and within each sector the dynamics drives to the maximally mixed state, so

    ρ(∞) = Σ_{w=0}^{N} p_w · (P_w / d_w),   p_w = Tr(P_w ρ₀),   d_w = C(N, w),

an N+1-way attractor, one basin per sector, each weighted by its untouched initial share.

The reading unifies two apparent exits of the dynamics: the lens regime (a single-excitation state, all weight in one sector) and the cusp regime (weight spread across several) are the same N+1-exit structure seen from different initial conditions. It is the excitation-number companion of the parity selection rule and of F4's stationary-mode count (kernel dimension N+1, which is F4's single-connected-component case of Π_c (|c| + 1), the same connectivity this theorem needs), and the proof runs in three independent steps: sector-population conservation (kinematic), the unique steady state per sector (dynamic), and the decay of off-diagonal coherence.

## Theorem

Let ρ₀ be any initial density matrix on N qubits. Let the dynamics be the Lindblad equation with Heisenberg Hamiltonian H = Σ J (X_k X_{k+1} + Y_k Y_{k+1} + Z_k Z_{k+1}) on any **connected** graph topology, and local Z-dephasing jump operators L_k = √(γ_k) Z_k with site-dependent rates γ_k > 0 at **every** site. Let P_w = Σ_{i: popcount(i)=w} |i⟩⟨i| be the projector onto the excitation-number sector w, and let d_w = C(N, w) be its dimension.

Then the asymptotic state is:

    ρ(∞) = Σ_{w=0}^{N}  p_w  ·  (P_w / d_w),        p_w = Tr(P_w ρ₀)

**In particular:** p_w(∞) = p_w(0) = Tr(P_w ρ₀). The asymptotic population of each excitation sector equals its initial population. The asymptotic state is a mixture of N+1 maximally-mixed sector states, weighted by the initial sector distribution.

---

## Proof

The theorem has two independent ingredients. Each is stated and proven separately.

### Step 1 (kinematic): sector population is a constant of motion

Define p_w(t) = Tr(P_w ρ(t)). The claim is dp_w/dt = 0 for all t, hence p_w(t) = p_w(0).

Let N_op = Σ_k (I − Z_k)/2 be the total-excitation-number operator. Direct computation:

1. **[H, N_op] = 0.** The Heisenberg bond XX + YY + ZZ commutes with N_op. ZZ is diagonal and trivially commutes. For XX+YY, note that XX+YY = 2(S_k⁺ S_{k+1}⁻ + S_k⁻ S_{k+1}⁺), which moves one excitation from site k to site k+1 or vice versa, preserving the total count.

2. **[L_k, N_op] = 0.** L_k = √(γ_k) Z_k is diagonal in the computational basis, as is N_op. Diagonal operators commute.

3. P_w is a spectral projector of N_op (the eigenspace projection 𝟙_{N_op = w}), so [H, P_w] = 0 and [L_k, P_w] = 0 as well.

4. The Lindblad generator applied to P_w from the left yields:
   dp_w/dt = Tr(P_w · L[ρ]) = −i Tr(P_w [H, ρ]) + Σ_k γ_k Tr(P_w (Z_k ρ Z_k − ρ))

   Using the cyclicity of the trace and [H, P_w] = 0:
   Tr(P_w [H, ρ]) = Tr([P_w, H] ρ) = 0.

   Using [Z_k, P_w] = 0 (diagonal commutation):
   Tr(P_w Z_k ρ Z_k) = Tr(Z_k P_w Z_k ρ) = Tr(P_w Z_k² ρ) = Tr(P_w ρ),
   so the dissipator contribution is also zero.

Therefore dp_w/dt = 0. This holds for all t, in particular p_w(∞) = p_w(0) = Tr(P_w ρ₀). End Step 1.

**Note:** Step 1 is equivalent to the sector conservation theorem proven independently in [Cusp-Lens Connection](../../experiments/CUSP_LENS_CONNECTION.md) and formalized in the [Parity Selection Rule](PROOF_PARITY_SELECTION_RULE.md). It is re-derived here for self-containedness.

Step 1 is also **wider than the theorem it feeds**. Nothing in it needs a rate to be nonzero or a graph to be connected: the computation above only uses that H and every L_k commute with N_op. So sector-population conservation holds for γ_k ≥ 0 on any graph. The stronger hypothesis enters at Step 2.

### Step 2 (dynamic): within each sector, the unique steady state is maximally mixed

The claim is: restricted to the (w, w) diagonal sector (operators of the form P_w X P_w), the Lindblad dynamics has exactly one steady state, namely P_w / d_w.

**Proof.** The fixed-point algebra of the restricted Lindblad generator L_w = −i[H_w, ·] + D_w is the intersection of (a) the decoherence-free subalgebra of D_w and (b) the commutant of H_w.

**(a) Decoherence-free subalgebra of D_w.** An operator A on the w-sector satisfies D_w[A] = 0 iff Z_k A Z_k = A for all k (with γ_k > 0). Since Z_k is diagonal with entries ±1, the condition Z_k A Z_k = A forces A[i,j] · (z_k(i) z_k(j) − 1) = 0 for all i, j, k where z_k(i) = 1 − 2·bit_k(i). For i ≠ j there exists at least one k where bit_k(i) ≠ bit_k(j) (since i ≠ j), giving z_k(i) z_k(j) = −1, hence A[i,j] = 0. Therefore the DFS consists of diagonal operators on the w-sector: a d_w-dimensional commutative algebra.

**(b) Commutant of H_w within diagonal operators.** A diagonal operator A satisfies [H_w, A] = 0 iff H_w[i,j] · (A[i,i] − A[j,j]) = 0 for all i, j. This requires A[i,i] = A[j,j] whenever H_w[i,j] ≠ 0. The XX+YY swap terms in the Heisenberg Hamiltonian connect any two w-excitation basis states that differ by moving one excitation to an adjacent site. On any connected graph, repeated swaps along edges reach every w-excitation configuration from any other: a connected graph has a spanning tree, and the transpositions along the edges of a tree on N vertices generate the full symmetric group S_N. Therefore all diagonal entries of A must be equal: A = α · I_w.

**(a) ∩ (b) = {α · I_w}.** The unique density matrix in this algebra is P_w / d_w. End Step 2.

**Numerical cross-check:** verified for N=3-7 chain, and at N=4-5 for the four connected topologies (chain, ring, star, complete) surveyed in [Symmetry Census](../../experiments/SYMMETRY_CENSUS.md). Zero exceptions.

### Step 2b (off-diagonal sectors decay to zero)

For the off-diagonal sector (w, w') with w ≠ w': an operator A = Σ c_{ij} |i⟩⟨j| with popcount(i) = w, popcount(j) = w' satisfies D[A] = 0 iff z_k(i) z_k(j) = 1 for all k with γ_k > 0, which, when every γ_k > 0, requires bit_k(i) = bit_k(j) for all k, hence i = j, hence popcount(i) = popcount(j), contradicting w ≠ w'. Therefore the DFS of the (w, w') sector is {0}. No operator survives dephasing, and the Hamiltonian cannot prevent this (it only redistributes, it cannot create fixed points that the dissipator destroys). Every off-diagonal sector decays to zero.

### Step 3 (assembly)

Combining Steps 1, 2, and 2b:

- All off-diagonal sector blocks of ρ(t) (coherences |w⟩⟨w'| with w ≠ w') decay to zero (Step 2b).

- Each diagonal sector block ρ_w(t) = P_w ρ(t) P_w converges to P_w/d_w scaled by its time-independent trace p_w(0) (Step 2, plus Step 1 which fixes the trace).

- Therefore ρ(∞) = Σ_w p_w(0) · (P_w/d_w), with p_w(0) = Tr(P_w ρ₀).

All three steps are now proven analytically. The theorem holds for any connected graph topology and any site-dependent Z-dephasing profile with all γ_k > 0.

End proof.

---

## Consequences

1. **Exit distribution is computable before evolution.** Given ρ₀, the asymptotic state is determined by N+1 numbers p_w = Tr(P_w ρ₀). No time evolution is required.

2. **Asymptotic state is a function of (p_0, ..., p_N) alone.** Two initial states with identical sector populations produce identical asymptotic states, regardless of their coherences or sector-internal structure. The vector (p_0, ..., p_N) is a complete invariant for the purpose of predicting ρ(∞).

3. **"Number of exits used" is the Hilbert-sector support size.** The number of attractors that ρ₀ contributes to is |{w : p_w > 0}|. A product state uses 1, a Bell+ state on two qubits uses 2, a GHZ state uses 2 (w=0 and w=N). The |+⟩⊗N state uses all N+1. Note: [Symmetry Census, "Which sectors are reachable"](../../experiments/SYMMETRY_CENSUS.md) counts these as Liouville blocks (w_bra, w_ket); a GHZ state populates 2 Hilbert sectors but 4 Liouville blocks (the two diagonal blocks plus the off-diagonal coherences that decay). Both counts are valid; they describe different objects.

4. **The cusp exit and lens exit are both instances of this theorem.** The lens exit (thermalization within a single excitation sector, see [Concentrator Geometry](../../experiments/CONCENTRATOR_GEOMETRY.md)) corresponds to initial states with p_1 = 1 (single-excitation sector only). The cusp exit (simultaneous thermalization across multiple sectors, see [Cusp-Lens Connection](../../experiments/CUSP_LENS_CONNECTION.md)) corresponds to initial states with p_w > 0 for two or more w. The "two exits" framing of earlier documents is a special case of the general N+1-exit structure.

---

## Scope and limitations

- **Holds for:** Heisenberg Hamiltonian on any **connected** graph topology, local Z-dephasing with γ_k > 0 at **every** site, any N ≥ 2. Step 1 alone is wider, holding for γ_k ≥ 0 on any graph.
- **Breaks for:** Amplitude damping (L_k = √(γ_k) σ_k⁻ does not commute with N_op), transverse-field Hamiltonians (H includes X_k or Y_k terms that do not commute with N_op), non-Markovian dynamics (no Lindblad form), and disconnected graphs (Step 2 part (b)'s adjacent-transposition argument has no path between two configurations in different components).
- **The argument needs a support covering every site; the conclusion need not.** With γ = 0 at some sites, Step 2 part (a) and Step 2b have no k to spend for pairs that agree wherever γ is nonzero, so the proof stops. The `scope` run of [The Blind Site](../../experiments/THE_BLIND_SITE.md) decides eight supports at N = 5 on the uniform chain, the full one and seven sparse ones, without an eigensolver. It bounds the undamped space of L, the span of every eigenvector whose eigenvalue lies on the imaginary axis (The Node Pair's peripheral space), from above by a GF(p) rank and from below by the N+1 sector projectors, which are stationary for every support; that is the instrument of [The Node Pair](PROOF_NODE_PAIR_RESOLVENT.md)'s gate G5, run on the full space. For seven of the eight supports the two bounds meet at 6: the undamped space is the span of the P_w, and every initial state converges to Σ_w p_w P_w/d_w. The eighth, the single centre seat, is where they do not meet at 6; there the upper bound is 24, and the next point meets it from below. Other N and other graphs are not decided here; the certificate runs on any support.
- **At the centre seat the long-time behaviour depends on where the state sits.** With M the site reversal, [Z_k, M] = 0 only at the fixed centre, so a centre support conserves the two mirror-parity projectors on top of the N+1 sector projectors. On the mirror-even side the same certificate, run on operators that live on the even space on both sides, meets at 6 again: a state with no weight in the mirror-odd space converges to the maximally mixed state on the even part of each sector. For |+⟩^5 that limit differs from Σ_w p_w P_w/d_w by exactly +1/48 on each mirror-symmetric configuration of the sectors w = 1–4, −1/192 on each configuration there that has a partner, and 5/192 on each partner coherence; the propagated state at t = 400 sits 3.962·10⁻¹⁰ from it in maximum entrywise difference, and 5/192 = 2.604·10⁻² from the per-sector state. The mirror-odd side need not mix at all. At popcount 1 it is spanned by (|0⟩−|4⟩)/√2 and (|1⟩−|3⟩)/√2 and has no amplitude at the centre, so the centre jump acts on it as the identity, while H acts as [[2,2],[2,0]]: a pure state there that is not an eigenvector of that matrix moves unitarily forever. That space is the centre's blind subspace in The Blind Site. The Node Pair §8 names, on the XY chain, the undamped space of the single-excitation block with light on the centre, End(D) ⊕ ⟨I_E⟩ (D the blind modes, I_E the identity on the others), of dimension m² + 1 at N = 2m + 1; its argument needs only that D is spanned by eigenvectors of H with a node at the centre, so the same space is undamped here, D being the mirror-odd space spanned above. The count closes at 24 = 6 + 10 + 8, each term in its own parity block, so the dimensions add: the 6 even parts of the sector projectors; 10 on the odd side, all of End(odd part) at popcounts 1 and 4 (4 each) and only the odd part's identity at popcounts 2 and 3; and 8 parity-crossing coherences, |00000⟩⟨o|, |o⟩⟨00000| with o in the popcount-1 odd part (The Blind Site §5, the (0,1) block of F152) and |11111⟩⟨o′|, |o′⟩⟨11111| with o′ in the popcount-4 odd part. The last eight agree on the centre, whose bit is 0 in every popcount-1 odd state and 1 in every popcount-4 one, and the `scope` run builds them as integer matrices and checks over ℚ, with no prime, that they have rank 8 and that their commutators with H stay in their span. That lower bound of 24 meets the GF(p) upper bound, so at the centre seat the undamped space is exactly this span, and the part of any state outside it decays.
- **Does not address:** Rates of convergence, transient structure, crossing behavior at CΨ = 1/4, or any dynamical feature on finite t. The theorem is asymptotic only.

---

## References

- [Symmetry Census](../../experiments/SYMMETRY_CENSUS.md) ("Asymptotic attractors per sector" and "Topology comparison": numerical verification of Steps 2–3)
- [Cusp-Lens Connection](../../experiments/CUSP_LENS_CONNECTION.md) (sector conservation theorem, precursor of Step 1)
- [Concentrator Geometry](../../experiments/CONCENTRATOR_GEOMETRY.md) (lens exit as p_1 = 1 case)
- [Parity Selection Rule](PROOF_PARITY_SELECTION_RULE.md) (formal proof of sector conservation, equivalent to Step 1)
- [The Blind Site](../../experiments/THE_BLIND_SITE.md) (what a dephasing support that misses sites leaves untouched, and the N = 5 case that bounds the scope above)

---

*Walked into existence by Tom and Claude (chat) on 2026-04-12 during Track B of the three-values investigation. Recorded so the experiment document ([Sector Projection Formula](../../experiments/SECTOR_PROJECTION_FORMULA.md)) could cite it rather than re-prove it.*
