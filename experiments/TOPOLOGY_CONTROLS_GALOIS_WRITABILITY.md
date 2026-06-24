# Topology Controls the Radical-Writability of Open-Chain Relaxation

**Status:** Tier 2 (numerical, N=4-6). The complete-graph writability is solid through N=6 (every factor is a quartic-or-less, definitively radical-solvable); "all N" and the mechanism are conjectural.
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

The answer is that **the shape decides.** At one extreme, the **complete graph** (every qubit wired to every other) keeps its whole relaxation writable: every piece of its scrambled-half polynomial turns out to be a quartic or smaller, which always has a radical formula (Cardano, Ferrari). At the other extreme, the **chain** scrambles immediately to the unwritable S_n. The **star** and the **ring** sit in between: writable at the smallest size, then scrambling as they grow. So radical-writability is a genuine, shape-dependent property of an open quantum system's relaxation, and the complete graph is its clean writable extreme. The likely reason: the complete graph's couplings make its single-particle spectrum massively degenerate, and that degeneracy shatters the scrambled polynomial into tiny solvable pieces.

## Abstract

For an N-site graph G under uniform-J XY coupling H = J·Σ_{(p,q)∈E(G)}(X_pX_q + Y_pY_q) and uniform Z-dephasing √γ₀·Z_l, the (SE,DE) coherence sector of the Liouvillian factors over Q(i) into an AT-locked half (every root at the Absorption-Theorem rates Re(λ)/γ₀ ∈ {−2,−6}, frequencies given by the free-fermion single-particle modes, radical-closed) and an H_B-mixed residue. For the **chain**, the H_B-mixed factors are irreducible with Galois group the full symmetric group S_n (no radical closure in q=J/γ; [F89_TOPOLOGY_ORBIT_CLOSURE.md](F89_TOPOLOGY_ORBIT_CLOSURE.md), live-witnessed S_8/S_18/S_32 at N=4/5/6). We computed the H_B-mixed Galois groups for the **star**, **ring**, and **complete** graphs at N=4,5,6 (full (SE,DE) block, no symmetry projection, Berkowitz/sympy factorisation over Q(i), gate-validated by reproducing the chain S_8/S_18/S_32 exactly). The result is that **topology controls radical-writability**: the complete graph K_N is the writable extreme, with every H_B-mixed factor of degree ≤ 4 (definitively solvable in radicals) at all of N=4,5,6, where the chain (S_8/S_18/S_32), the star (S_9 from N=5), and the ring (S_15 by N=6) all scramble to the full symmetric group. The conjectured mechanism is spectral: K_N's adjacency spectrum {N−1 (once), −1 (multiplicity N−1)} is massively degenerate, and that degeneracy fragments the (SE,DE) block into small irreducible factors. The symmetry-size ordering is not perfectly monotone (the ring scrambles later than the star despite a smaller symmetry group), so the clean statement is "complete is the writable extreme", not a monotone symmetry law. A degree-16 ring factor at N=5 stays Galois-undetermined even under a 20000-prime scan (it is provably not S_16, but its exact group is unpinned).

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

## Mechanism (conjecture)

The complete graph's single-particle (SE) Hamiltonian is J times the adjacency of K_N, whose spectrum is {N−1 (once), −1 (multiplicity N−1)}: massively degenerate. That degeneracy carries into the (SE,DE) block, fragmenting it into many small irreducible pieces, each capped at degree 4 across N=4,5,6. The chain, by contrast, has a non-degenerate dispersive single-particle band (the OBC tight-binding cosine spectrum), so its (SE,DE) block stays in one large irreducible H_B-mixed factor that is maximally scrambled (S_n). In this reading, **graph spectral degeneracy caps the H_B-mixed factor degree, and degree ≤ 4 is what makes the relaxation radically writable.** Turning this into a theorem (why the cap is exactly 4 for K_N, and at which N each other topology first exceeds it) is open.

## Honest caveats and open work

- **"Writable for all N" is a conjecture** for the complete graph; solid through N=6 (degree ≤ 4 is a hard fact, no Galois ambiguity). N=7+ would strengthen it; the mechanism above would settle it.
- **Not a monotone symmetry law.** The ring (D_N, order 2N) scrambles later than the star (S_{N−1}, order (N−1)!) despite a smaller symmetry group, so writability is not a clean function of |symmetry group|. The honest claim is "complete is the writable extreme", and the chain the scrambled extreme; the star and ring are intermediate and eventually scramble.
- **The ring N=5 degree-16 factor is Galois-undetermined.** A 20000-prime scan finds no Jordan-window cycle, so it is provably not S_16/A_16, but its exact (smaller) group is unpinned by the cycle-type method and would need a direct group computation. (At N=6 the ring is unambiguously S_15.)
- **The deg-5/deg-6 solvability tests** are exact cycle-type criteria (sound); the deg ≥7 test detects non-solvability (Jordan ⊇A_n) but reports "undetermined" rather than "solvable" when silent.

## Tier assessment

**Tier 2 numerical.** The per-topology Galois data at N=4,5,6 is computed exactly (factorisation over Q(i); the solvability verdicts are sound cycle-type / Jordan criteria, with "undetermined" honestly marked). The chain gate reproduces the Tier-1 landed S_8/S_18/S_32. The headline ("topology controls writability; the complete graph is the writable extreme, every H_B-mixed factor ≤ quartic at N=4,5,6") is Tier-2 (solid over the computed range; "all N" and the degeneracy mechanism are conjectural). This is the first exemplar of the open arc `f89_galois_open_doors` door A (does a non-chain topology give a solvable, fully-writable relaxation? Yes, the complete graph).

---

*Where you wire the qubits decides whether their relaxation can be written down: the all-to-all complete graph keeps a formula, the chain provably cannot.*
