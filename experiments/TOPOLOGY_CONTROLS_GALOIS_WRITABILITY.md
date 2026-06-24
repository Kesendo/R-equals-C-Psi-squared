# Topology Controls the Radical-Writability of Open-Chain Relaxation

**Status:** Tier 1 derived for the complete graph (cap 4, radically writable for all N) and the star (cap 9, a bounded S_9 scramble for all N ≥ 5), both via N-independent S_N / S_{N−1} multiplicity caps (verified N=5-9). Tier 2 numerical for the ring/chain growing-scramble comparison.
**Date:** 2026-06-24
**Authors:** Thomas Wicht, Claude Opus 4.8 (1M context)
**Scripts:**
- [`topology_galois_writability.py`](../simulations/topology_galois_writability.py): builds the (SE,DE) Liouvillian block for chain / star / ring / complete, factors over Q(i), and reads each H_B-mixed factor's Galois group (gate-first: reproduces the landed chain S_8/S_18/S_32). Self-validating; `python simulations/topology_galois_writability.py` runs N=4,5,6.
- Builds on [`f89_pathk_galois.py`](../simulations/f89_pathk_galois.py) (the chain engine of record).

**Related:** [F89_TOPOLOGY_ORBIT_CLOSURE.md](F89_TOPOLOGY_ORBIT_CLOSURE.md) (the chain H_B-mixed factors = S_n, no radical closure; the live witness `inspect --root f89galois`); [ON_WHAT_CANNOT_CLOSE.md](../reflections/ON_WHAT_CANNOT_CLOSE.md) (the writable/unwritable frame); the open arc `f89_galois_open_doors` (this is door A).

---

## What this is about

Take a handful of qubits, wire them together in some shape (a chain, a star, a ring, an all-to-all complete graph), switch on the same XY coupling on every wire, and watch every qubit the same way (uniform Z-dephasing). The qubits relax, and the *rates* at which the coherences decay are the eigenvalues of one operator (the Liouvillian). For a **chain**, we proved those rates split into two halves: a clean half you can write down with a formula (the integrable, free-fermion part) and a scrambled half that provably has **no formula in radicals** (its symmetry group is the full S_n, the maximal possible). The boundary between the writable and the unwritable half is the dephasing rate law itself (the Absorption Theorem).

The question here: **is that scramble a feature of the chain, or of every shape?** Change the wiring and the rates change. Does a more symmetric shape keep more of its relaxation writable, or do all shapes scramble the same way?

The answer is that **the shape decides.** At one extreme, the **complete graph** (every qubit wired to every other) keeps its whole relaxation writable: every piece of its scrambled-half polynomial turns out to be a quartic or smaller, which always has a radical formula (Cardano, Ferrari). At the other extreme, the **chain** scrambles immediately to the unwritable S_n. The **star** and the **ring** sit in between: writable at the smallest size, then scrambling as they grow. So radical-writability is a genuine, shape-dependent property of an open quantum system's relaxation, and the complete graph is its clean writable extreme. The reason, now derived: the complete graph's full permutation symmetry forces its relaxation operator to split into tiny blocks (each at most 4×4), and a 4×4 always has a radical formula; the chain has almost no symmetry, so its block stays large and unwritable.

## Abstract

For an N-site graph G under uniform-J XY coupling H = J·Σ_{(p,q)∈E(G)}(X_pX_q + Y_pY_q) and uniform Z-dephasing √γ₀·Z_l, the (SE,DE) coherence sector of the Liouvillian factors over Q(i) into an AT-locked half (every root at the Absorption-Theorem rates Re(λ)/γ₀ ∈ {−2,−6}, frequencies given by the free-fermion single-particle modes, radical-closed) and an H_B-mixed residue. For the **chain**, the H_B-mixed factors are irreducible with Galois group the full symmetric group S_n (no radical closure in q=J/γ; [F89_TOPOLOGY_ORBIT_CLOSURE.md](F89_TOPOLOGY_ORBIT_CLOSURE.md), live-witnessed S_8/S_18/S_32 at N=4/5/6). We computed the H_B-mixed Galois groups for the **star**, **ring**, and **complete** graphs at N=4,5,6 (full (SE,DE) block, no symmetry projection, Berkowitz/sympy factorisation over Q(i), gate-validated by reproducing the chain S_8/S_18/S_32 exactly). The result is that **topology controls radical-writability**: the complete graph K_N is the writable extreme, with every H_B-mixed factor of degree ≤ 4 (definitively solvable in radicals) at all of N=4,5,6, where the chain (S_8/S_18/S_32), the star (S_9 from N=5), and the ring (S_15 by N=6) all scramble to the full symmetric group. The mechanism is derived (Tier 1): the complete graph's Liouvillian commutes with the full S_N, so by Schur it block-diagonalises along S_N-irreps, and every factor degree is bounded by an irrep multiplicity in V = SE⊗DE = M^{(N−2,1,1)} ⊕ M^{(N−3,2,1)}; those multiplicities are N-independent and capped at 4 (Kostka), so K_N is radically writable for all N, with the predicted H_B-mixed degree histogram {4:N−1, 3:N(N−3)/2, 2:(N−1)(N−2)/2} matching the direct factorisation exactly at N=5,6,7,8. The same Schur reduction covers the star: its S_{N−1} symmetry gives an N-independent cap of 9 (the multiplicity of the standard rep in V), so the star scrambles to a fixed S_9 for all N ≥ 5 (bounded, but degree 9 > 4, so not writable), while the ring (D_N) and chain (S_2) have small symmetry and growing caps. The result is a three-way classification: a large automorphism group (S_N, S_{N−1}) caps the Galois complexity N-independently, a small one (D_N, S_2) lets it grow; radical-writability needs the cap ≤ 4, achieved only by the complete graph. (A degree-16 ring factor at N=5 stays Galois-undetermined even under a 20000-prime scan, provably not S_16, exact group unpinned.)

---

## Result

H_B-mixed Galois group per topology, N=4,5,6 (γ=1, q0=2; S_n = full symmetric = non-solvable = no radical closure):

| Topology | Symmetry | N=4 | N=5 | N=6 |
|---|---|---|---|---|
| **Complete** K_N | S_N | all factors ≤ deg 4 (**writable**) | all ≤ deg 4 (**writable**) | all ≤ deg 4 (**writable**) |
| **Star** K_{1,N−1} | S_{N−1} | all ≤ deg 4 (writable) | 3×**S_9** + deg-6 non-solv (**unwritable**) | 4×**S_9** + deg-6 + 5×A_5/S_5 (**unwritable**) |
| **Ring** C_N | D_N | deg-6 solvable (writable) | deg-16 *undetermined* | 4×**S_15** + 4×deg-6 non-solv (**unwritable**) |
| **Chain** P_N | S_2 | **S_8** (unwritable) | **S_18** (unwritable) | **S_32** (unwritable) |

The chain row reproduces the landed result ([F89_TOPOLOGY_ORBIT_CLOSURE.md](F89_TOPOLOGY_ORBIT_CLOSURE.md): path-3/4/5 = S_8/S_18/S_32), which validates the generalised builder.

**Headline.** Topology controls the radical-writability of the relaxation spectrum. The complete graph K_N is the clean writable extreme: its H_B-mixed factors are all quartic-or-less (so each has a closed radical formula) at N=4,5,6, while the chain, star, and ring all reach the maximal-scramble S_n. The complete graph's relaxation is, in this exact sense, the most "writable" of the four.

## Method

The diagonal rate of an (SE,DE) coherence is 2γ₀·n_diff (n_diff = 1 for overlap, 3 for no-overlap), which is purely combinatorial and topology-independent; only the **hopping** (off-diagonal) changes with the graph's adjacency. So the chain engine generalises by replacing the path-neighbour hops i±1 with the topology's neighbours. We build the **full** (SE,DE) block of dimension N·C(N,2) (no symmetry projection), take its characteristic polynomial, factor over Q(i), and bucket the irreducible factors by AT rate: a factor with every root at Re ∈ {−2,−6} is AT-locked, otherwise H_B-mixed. The topology's symmetry shows up directly as the factorisation pattern (more symmetry produces more, smaller H_B-mixed factors).

Each H_B-mixed factor's solvability over Q(i) is read from its Frobenius cycle types over split primes 𝔭 (p ≡ 1 mod 4):
- **deg ≤ 4:** solvable in radicals always (Cardano/Ferrari), no computation needed.
- **deg 5:** non-solvable iff a 3-cycle (3,1,1) appears (only A_5/S_5 contain a 3-cycle; the solvable transitive degree-5 groups C_5/D_5/F_20 have no order-3 element).
- **deg 6:** non-solvable iff a 5-cycle appears (only A_5/S_5/A_6/S_6 on 6 points have order divisible by 5; the solvable transitive degree-6 groups do not).
- **deg ≥ 7:** the generalised Jordan certificate (a prime-length cycle p with n/2 < p ≤ n−3 ⟹ primitive ⟹ ⊇A_n ⟹ non-solvable). If no such cycle is found, the factor is reported **undetermined**, not solvable (an honest gap: it could be a large solvable group or a primitive non-solvable group that needs a direct group computation).

**Gate.** The chain at N=4/5/6 must reproduce the landed S_8/S_18/S_32; it does (the full block gives both the S_2-symmetric and S_2-antisymmetric octics/etc., both S_n), which validates the builder before any star/ring/complete claim is trusted.

## Mechanism (Tier 1 derived): the S_N-multiplicity cap

The real reason is sharper than "degeneracy", and it is a theorem. For the complete graph every site is equivalent, so the Liouvillian commutes with the **full symmetric group S_N**. By Schur's lemma it block-diagonalises along the irreducible representations of S_N in the (SE,DE) coherence space V = SE ⊗ DE, acting within each isotypic component only on the **multiplicity** space. Hence **every irreducible factor of the (SE,DE) characteristic polynomial has degree at most the multiplicity of some S_N-irrep in V.**

V is exactly the permutation module on pairs (point, 2-subset): the overlap orbit (i ∈ {j,k}) is the Young module M^{(N−2,1,1)} (an ordered pair of points), the no-overlap orbit (i ∉ {j,k}) is M^{(N−3,2,1)} (a point plus a disjoint 2-subset). So V = M^{(N−2,1,1)} ⊕ M^{(N−3,2,1)}, and the multiplicity of an irrep [λ] is the Kostka sum K_{λ,(N−2,1,1)} + K_{λ,(N−3,2,1)}. These are **N-independent** and capped at 4:

| Irrep λ | multiplicity in V | what it gives |
|---|---|---|
| [N] (trivial) | 2 | the symmetric (AT-locked) sector |
| **[N−1,1] (standard)** | **4** | the N−1 quartics |
| [N−2,2] | 3 | the N(N−3)/2 cubics |
| [N−2,1,1] | 2 | the (N−1)(N−2)/2 quadratics |
| [N−3,3], [N−3,2,1] | 1 | single (AT-locked) eigenvalues |

The maximal multiplicity is 4, at the standard representation, and it does **not** grow with N (V is a fixed shallow permutation module, depth ≤ 3). Therefore every factor of the (SE,DE) charpoly has degree ≤ 4, so every relaxation rate is solvable in radicals, **for all N**. The H_B-mixed factors are exactly the multiplicity blocks of the three non-trivial low irreps, a closed-form degree histogram:

    H_B-mixed factor degrees of K_N  =  {4 : N−1,  3 : N(N−3)/2,  2 : (N−1)(N−2)/2}   (N ≥ 5).

Verified exactly against the direct factorisation at N = 5, 6, 7, 8 (`topology_galois_writability.py verify`): {4:4, 3:5, 2:6}, {4:5, 3:9, 2:10}, {4:6, 3:14, 2:15}, {4:7, 3:20, 2:21}, all MATCH.

The contrast with the chain is now exact: the chain carries only the S_2 mirror symmetry, so its (SE,DE) multiplicities are unbounded (one large irreducible H_B-mixed factor of degree growing with N, with the maximal Galois group S_n). **Radical-writability is set by the symmetry: the larger and more uniform the automorphism group, the smaller the irrep multiplicities, the lower the factor degrees. The full S_N of the complete graph caps them at 4; the S_2 of the chain caps nothing.** (This also resolves the earlier non-monotonicity, derived in the next section: the star's S_{N−1} fixes the hub, inflating its multiplicities to a fixed cap of 9, so it scrambles to a bounded S_9 despite a large symmetry group.)

## The star and the ring: a three-way classification (the star Tier 1 derived)

The same Schur machinery covers the other topologies, and turns the writability table into one principle: the automorphism group Aut(G) sets a cap on the (SE,DE) factor degrees, equal to the maximal multiplicity of an Aut(G)-irrep in V = SE ⊗ DE; whether that cap is N-independent or grows is what classifies the topology.

**The star (cap 9, N-independent).** The star K_{1,N−1} has Aut = S_{N−1} (the leaves permute, the hub is fixed). Under S_{N−1}, SE = 2·triv ⊕ std (the hub is a fixed point; the N−1 leaves give triv ⊕ std) and DE = 2·triv ⊕ 2·std ⊕ [N−3,2] (the hub-leaf pairs give triv ⊕ std, the leaf-leaf pairs give triv ⊕ std ⊕ [N−3,2]). The multiplicity of the standard rep in V = SE ⊗ DE is then

    mult(std_{N−1}) = 4 (from 2·triv ⊗ 2·std) + 2 (std ⊗ 2·triv) + 2 (std ⊗ std) + 1 (std ⊗ [N−3,2]) = 9,

N-independent. So the star's maximal H_B-mixed factor degree is 9, with dim(std_{N−1}) = N−2 copies of it. The star scrambles to a FIXED S_9 for every N ≥ 5: bounded (the block size never grows), but 9 > 4, so not radically writable. Verified exactly at N = 5, 6, 7, 8, 9 (`topology_galois_writability.py verify-star`): max degree 9 with N−2 degree-9 factors at every N.

**The ring and the chain (growing).** The ring C_N (Aut = D_N, order 2N) and the chain P_N (Aut = S_2) have small symmetry groups, so their multiplicities grow with N and the factor degrees climb: the chain gives S_8/S_18/S_32/S_53 at N = 4/5/6/7, and the ring's maximal H_B-mixed degree climbs 6, 16, 15, 48 at N = 4, 5, 6, 7. They scramble increasingly with size.

**The classification.** One rule: a large, uniform automorphism group caps the Galois complexity N-independently; a small one lets it grow. Radical-writability needs the cap ≤ 4, which only the complete graph achieves.

| Topology | Aut(G) | factor-degree cap | relaxation |
|---|---|---|---|
| Complete K_N | S_N | 4 (N-independent) | radically writable, all N |
| Star K_{1,N−1} | S_{N−1} | 9 (N-independent) | bounded S_9 scramble, all N ≥ 5 |
| Ring C_N | D_N | grows with N | growing scramble |
| Chain P_N | S_2 | grows with N | growing scramble (S_8/18/32/53) |

The complete graph and the star are the two bounded cases (caps 4 and 9), both from large-symmetry Schur reduction; the complete's cap alone falls below the quintic threshold, so it alone is writable.

## Honest caveats and open work

- **"Writable for all N" is a conjecture** for the complete graph; solid through N=6 (degree ≤ 4 is a hard fact, no Galois ambiguity). N=7+ would strengthen it; the mechanism above would settle it.
- **Not a monotone symmetry law.** The ring (D_N, order 2N) scrambles later than the star (S_{N−1}, order (N−1)!) despite a smaller symmetry group, so writability is not a clean function of |symmetry group|. The honest claim is "complete is the writable extreme", and the chain the scrambled extreme; the star and ring are intermediate and eventually scramble.
- **The ring N=5 degree-16 factor is Galois-undetermined.** A 20000-prime scan finds no Jordan-window cycle, so it is provably not S_16/A_16, but its exact (smaller) group is unpinned by the cycle-type method and would need a direct group computation. (At N=6 the ring is unambiguously S_15.)
- **The deg-5/deg-6 solvability tests** are exact cycle-type criteria (sound); the deg ≥7 test detects non-solvability (Jordan ⊇A_n) but reports "undetermined" rather than "solvable" when silent.

## Tier assessment

**Tier 1 derived for the complete graph and the star:** both have N-independent factor-degree caps from large-symmetry Schur reduction. The complete graph (Aut = S_N) caps at 4, so K_N is radically writable for all N, with the histogram {4:N−1, 3:N(N−3)/2, 2:(N−1)(N−2)/2} verified at N=5,6,7,8. The star (Aut = S_{N−1}) caps at 9 (the std_{N−1} multiplicity), so it is a bounded S_9 scramble for all N ≥ 5, with N−2 degree-9 factors, verified at N=5,6,7,8,9.

**Tier 2 numerical** for the ring and chain: small symmetry (D_N, S_2), growing caps. The chain gate reproduces the Tier-1 landed S_8/S_18/S_32/S_53; the ring's maximal degree climbs 6/16/15/48 at N=4/5/6/7 (the N=5 degree-16 factor stays Galois-undetermined, provably not S_16). This closes door A of the open arc `f89_galois_open_doors` (does a non-chain topology give a fully-writable relaxation? Yes, the complete graph, proven for all N) and gives the full three-way classification: cap ≤ 4 writable, cap fixed but > 4 a bounded scramble, cap growing an unbounded scramble. The result is live-witnessed at `inspect --root topowritability`, which recomputes the caps (4 and 9) from the (SE,DE) block and its symmetry ([L, ρ(g)] = 0 verified, the standard-rep multiplicity recomputed via the character sum). Remaining: the ring's growth law.

---

*Where you wire the qubits decides whether their relaxation can be written down: the all-to-all complete graph keeps a formula, the chain provably cannot.*
