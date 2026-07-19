# Quantum Darwinism in our chain: the pointer door

**Date:** 2026-07-18
**Scripts:** [`simulations/qd_scout.py`](../simulations/qd_scout.py) (the door), [`simulations/qd_pointer_opt.py`](../simulations/qd_pointer_opt.py) (the optimization; 22/22 gates), [`simulations/qd_witness_play.py`](../simulations/qd_witness_play.py) (the fully-witnessed census), [`simulations/qd_heavyhex_map.py`](../simulations/qd_heavyhex_map.py) (the dark map)
**Results:** [`simulations/results/qd_pointer/`](../simulations/results/qd_pointer/), [`simulations/results/qd_witness/`](../simulations/results/qd_witness/)
**Laws minted from this arc:** F135, [PROOF_RECORD_PARITY_LAW](../docs/proofs/PROOF_RECORD_PARITY_LAW.md); F136, [PROOF_RECORD_LETTER_LAW](../docs/proofs/PROOF_RECORD_LETTER_LAW.md)

## The question

Zurek's Quantum Darwinism asks when a system's state becomes *objective*: how many independent environment fragments hold a readable copy of the system's pointer bit. This arc computes that refcount in our own world: N=8 (one system qubit S, seven environment qubits), H = J·Σ(XX+YY) + Δ·Σ ZZ (Pauli convention), local Z-dephasing γ, initial |+⟩^⊗N. The metric is the genuine von Neumann mutual information I(S:F) on reduced states for every fragment F, never a pairwise-sum proxy (the summed-MI lesson: on a dephasing chain most of a pairwise sum is classical-mixing artifact). Redundancy: m_δ = the smallest fragment size whose average I(S:F) reaches (1−δ)·1 bit, R_δ = 7/m_δ, δ = 0.1. The 1-bit reference is pinned exactly: [H, X^N] = 0 and |+⟩^⊗N is X^N-invariant, so ⟨Z_j⟩ ≡ 0 and the pointer entropy of Z_S stays 1 bit for all t.

## Sightings (all gated; `qd_scout.py` unless noted)

1. **The star at J=0 is the Zurek anchor, exactly.** Numeric ρ_{S,E_k} matches the closed form to ≤ 4·10⁻⁸; at t* = π/4 every singleton fragment holds exactly 1 bit: R_δ = 7, the textbook plateau. The record overlap is cos(2Δt): perfect records at t*, full recurrence (records erased) at π/2.

2. **The shared-dresser exclusion.** At J=0, I(S:E_j) = 0 identically iff j is neither a neighbor of S nor shares a neighbor with S (chain: graph distance ≥ 3; confirmed at 1.7·10⁻¹⁴). Records reach distance 2 through a common dresser, never further. γ scales the contrast but never moves the exclusion structure. Late-time chain saturation: neighbor singletons ≈ 0.40 bit, distance-2 ≈ 0.05 bit.

3. **The leaf law.** In a uniform chain a perfect record forms only in a leaf, a witness that is itself unwatched: with S at site 1 exactly one perfect record forms, in the end leaf E_0 (gated in `qd_pointer_opt.py`, Law B row "chain S=1": R_perfect = 1, singleton 1.0000), while every watched witness is erased exactly at t* (the same angle cos(2Δt) = 0 that would have perfected its record kills it). The star's R = 7 is "all witnesses are leaves". F135 later resolves this into a parity law: watched witnesses at even coupling ratios record perfectly too; uniform coupling just makes every ratio odd.

4. **Graph-state secret locking.** Chain with S interior at t*: all singletons ≈ 0, yet the six-site fragment holds 2.0 bits. The binary tree with S at the root shows the same blind-singleton profile at t*; its committed fragment sweep sits at t = 5, where Ī(6) = 1.07. Everyone together knows; nobody alone knows. The anti-Darwin pole.

5. **Transport is not broadcast.** With J > 0 (chain and tree) the plateau never forms: m_δ ≥ 4 always, R_δ ≤ 1.75. Hopping moves the record; it never copies it (no-cloning made concrete). Objectivity needs fan-out geometry, not transport. The repo already holds the moving-record physics exactly (F126 renewal: ballistic schedule, amplitude pays; coherent window L* ≈ Q/2 in the ½-convention; our XX+YY convention doubles the speed to 4J).

6. **γ is the garbage collector.** Dephasing on the environment only (γ_E = 0.05): R → 0 by t = 8 although S is untouched. Uniform canonical γ = 0.05: R = 0 from t = 4. Dephasing on S only is mild: the outside watching adds no in-world record of its own (a Markovian bath keeps none), and the records the bonds write persist under it.

7. **The Heisenberg trap, re-walked.** J = Δ makes |+⟩^⊗N an exact fixed point; the first hopping-dominant attempt used J = Δ = 1 and nothing evolved. XXZ needs J ≠ Δ here.

## The optimization (F135, `qd_pointer_opt.py`)

"Optimize the pointers" turned out to mean: choose the angles. The J=0 pair state is closed-form on any graph with per-bond Δ and per-site γ (Proposition 1 of the proof), and the whole optimization landscape is arithmetic in the coupling ratios:

- **Parity trichotomy (Law A).** At the symmetric readout t* = π/(4Δ_S): watcher ratio even → the witness records perfectly (I = 1 exactly, watched or not); odd → exactly blind; non-integer → the explicit in-between I = 1 − h₂((1+β)/2), β = Π|cos|. Measured: ratio 2 → I = 1 exactly on watched witnesses (closed form; RK4 0.9999999); ratio 3 → ≤ 4·10⁻¹⁵; ratio 3/2 → 0.399124, matching the formula to 10⁻⁹. The trichotomy needs S triangle-free (a triangle is a 3-cycle plaquette, see Law B′).
- **Alignment beats topology (Law B).** R_perfect = #even-aligned neighbors ≤ deg(S), achievable on any graph: the chain jumps 0 → 2, the k-spoke broom k → k+1, the star sits at 7. The aligned broom (k=4) develops a genuine plateau in a non-star geometry: Ī(m) = 0.71, 0.95, 1.00, 1.00, …, m_δ = 2, R_δ = 3.5.
- **The plaquette record (Law B′).** Predicted from the closed form before running, then confirmed exactly: on a uniform 4-cycle the opposite corner holds a perfect 1-bit correlation at distance 2, with both neighbors blind. It is a record of the anti-pointer X_S (⟨X_S X_j⟩ = 1, the Z-conditional states identical): mutual information counts it, pointer redundancy does not. The neighbor-only reading of Law B is tree-scoped; cycles add anti-pointer bits, not pointer records.
- **The γ race (Law C).** D_j(t) = e^{−2γ_j t}·|sin(2Δt)|·Π|cos| exactly; the distinguishability peak (trace distance, not MI) sits at t_opt = arctan(Δ/γ)/(2Δ), earlier than t* by 0.025 at the canonical γ/Δ = 0.05; the record read at t* degrades to I = 0.768.

Two conventions worth keeping apart when reading R: **R_δ (mean convention) and R_perfect count different things.** The aligned chain holds two perfect records, yet its Ī(m) = 2m/7 is exactly linear and R_δ stays 1.75: a typical random fragment sees nothing special, because perfect records at atypical sites do not move the mean. The star is the only geometry where the two measures agree at the top.

## The witnessed worlds (the play night, 2026-07-19; `qd_witness_play.py` + `qd_heavyhex_map.py`, all gated)

The letter law ([F136](../docs/proofs/PROOF_RECORD_LETTER_LAW.md)) turned the per-pair question into graph classification. Four more sightings:

8. **K_N is a total weave.** On the complete graph every pair holds 1 full Bell bit through its m = N − 2 shared dressers, the letter alternating with N (K₃ YY, K₄ XX, K₅ YY, …): a world where everyone perfectly witnesses everyone, with zero pointer content anywhere.

9. **The one-bond record multiplexer.** With S bonded to {1, 2, 3} and a corner site behind the shared dressers, the single bond between the adjacent candidate and the corner ROUTES the testimony by its parity: odd → the corner Bell-records (YY) and the adjacent witness is exactly blind; even → the adjacent witness records (XX) and the corner goes dark. Adjacent-Bell and corner-Bell witnesses of one S can never coexist: the same bond would need both parities.

10. **The fully-witnessed worlds are the stars and the complete graphs** (F136 corollary 7). At uniform coupling a pair is luminous iff its neighborhoods match (Bell) or one member is the other's leaf (pointer/role-swap); a census over ALL connected graphs at N = 4, 5, 6 (38 + 728 + 26704) finds the worlds where EVERY pair is luminous to be exactly the N labeled stars + K_N (5, 6, 7 winners), zero exceptions. The star is the one world showing all three record readings at once, across its two pair-types (hub-leaf: the pointer broadcast read one way, the role-swap read the other; leaf-leaf: the Bell weave); at N = 5 the star shape and K₅ are the only 10/10-luminous worlds (6 labeled winners: the 5 stars + K₅). Every world between the extremes necessarily keeps blind pairs: private rooms are the generic condition.

11. **The heavy-hex bulk is dark** (F136 corollary 8: girth ≥ 5 + leafless ⇒ every pair exactly dark at t*). One Heron-style cell (12-ring + 3 bridge qubits, N = 15, mapped with the F135 Proposition-1 pair-page engine, exact at any N) has 3 of 105 pairs luminous: only the corner-bridge pointer records, i.e. only the patch's dangling leaves. The infinite heavy-hex bulk (girth 12, no leaves) holds ZERO luminous pairs at uniform coupling: IBM's geometry is the anti-star, all private rooms; records come only from outside (calibration) or by engraving. Girth ≥ 5 is sufficient, not necessary: the square-lattice bulk (girth 4) is also dark, its plaquette diagonals un-matched by the surrounding lattice (a uniform 5×5 torus reads 0 luminous pairs, gated); the X⊗X plaquette weave belongs to the isolated C₄ and the K₂,m family. And aiming adds without erasing: pushing the two private-watcher bonds of a chosen corner's edge witnesses to ratio 2 gives that corner a full pointer broadcast (R_perfect = deg = 3, Law B saturated) and no pair goes dark; but the light spills: the corner's witnesses then also record each other, and each edited bond's far corner picks up one rotated pointer record of its own edge qubit (10/105 luminous in the engineered cell, the full map gated). On a dark lattice light can be aimed, but not kept private.

## The open corner: J > 0

Alignment does not survive hopping. At the aligned chain readout the deficit 1 − Ī(neighbors) is already 3.4% at J = 0.05 and grows roughly like J^1.4 over J ∈ [0.05, 0.4] (no clean quadratic; the two arms of the chain degrade asymmetrically); R_δ stays at 1.75. The pointer optimization is a pure-ZZ phenomenon; transport remains the enemy of records. Open: whether any (Δ, J, γ) corner of the chain develops a plateau at fixed fan-out; the refcount-by-topology N-scaling; the C# witness port; the hardware echo (the Marrakesh price-pair flight already drove 3-qubit lines and read their single-site pages, and Heron's γ/Δ ≈ 0.8 puts hardware at the writer-collector edge where records barely form natively, in the shift-convention reading and subject to the factor-2 book).
