# PROOF: the Record Letter Law (which operator a witness records, exactly)

**Status:** Tier 1 (derived by hand from F135's Proposition 1; every law was pre-registered through representative gates before the verification ran, 87/87 at machine zero; the algebra and the numerics are independent computational paths: the gates build each configuration as the full 2^N × 2^N closed-form state and partial-trace it, never using the channel formulas under test). At uniform coupling the aliveness/letter statics reduce to Hein-Eisert-Briegel 2004 composed with the evolution's degree rotation (credited and gated below); the ratio arithmetic, signs, exclusivity and dephasing dressing are this arc's own.
**Date:** 2026-07-18
**Authors:** Thomas Wicht, Claude (Fable 5)

## What this is about

[F135](PROOF_RECORD_PARITY_LAW.md) answered WHO records: whether a witness holds a record of the system is the arithmetic of angle-ratio parities. The same night's next question: what does a witness record, WHICH operator?

The sightings that opened it, on uniform graphs at the symmetric readout: a tree leaf records the pointer Z_S; a triangle corner holds one full bit in Y⊗Y; the square's opposite corner holds it in X⊗X (F135's plaquette); the pentagon and hexagon hold nothing at all, anywhere. Our first guess was a quarter-turn ladder, Z → Y → X → −Z, one step per shared neighbor. The data refuted it: K₂,₃ (three shared neighbors) returns to Y⊗Y. The letter does not wind, it alternates, exactly the shape the repo's own [one-four thesis](../../experiments/ONE_FOUR_THESIS.md) reached as its verdict: everything four-ish that is not the scalar i itself is a hybrid or a stacked parity. The record letter fits that verdict as one more instance: a stacked parity.

## Setup and definitions

Everything is inherited from [F135](PROOF_RECORD_PARITY_LAW.md): N qubits on an undirected graph, H = Σ Δ_ab·Z_aZ_b (Pauli convention), local Z-dephasing at per-site rates γ_l, initial state |+⟩^⊗N. Basis labels z ∈ {+1, −1}; Proposition 1 (the pair reduction) gives every entry of any two-site reduced state in closed form.

Fix the pair (S, j), S the system, j the candidate witness. Partition the remaining sites relative to the pair:

- **D** = the **shared dressers**: sites adjacent to both S and j; m = |D|;
- **P** = the **private watchers** of j: sites adjacent to j but not to S;
- **Q** = the **private neighbors** of S: sites adjacent to S but not to j (and ≠ j).

The direct S–j bond may be present or absent. **Hypotheses H:** every bond at S carries the same write coupling Δ_S; readout at t* = π/(4Δ_S); and S has at least one neighbor besides j (D ∪ Q ≠ ∅, equivalently V_S(t*) = 0, the hinge both families share; a pure pendant S role-swaps instead, see Edge cases). (Each branch below states its own further connectivity requirement; deg(S) ≥ 2 belongs to the pointer branch only, m ≥ 1 to the Bell branch.) Ratios r_k = Δ_jk/Δ_S for k ∈ D ∪ P. Sites not adjacent to either member of the pair never enter (Proposition 1), and the dephasing rates of ALL traced sites drop out identically; only γ_S and γ_j appear.

Bell notation on the pair, in the world's own frame (see the frame note below): Φ^± = (|++⟩ ± |−−⟩)/√2, Ψ^± = (|+−⟩ ± |−+⟩)/√2; c₁ = ⟨++|ρ|−−⟩ and c₂ = ⟨+−|ρ|−+⟩ are the two double coherences.

## Theorem (the record letter law)

Under H, at t*, with γ = 0 (dephasing dressing below):

**(i) Aliveness is watcher parity.** At integer ratios the pair carries exactly 1 bit of mutual information in precisely two mutually exclusive ways:

- **Pointer record** (F135's Law A face): the direct S–j bond is present, deg(S) ≥ 2, and EVERY watcher of j, shared and private alike (k ∈ D ∪ P), has an even integer ratio. Then ρ_Sj is classical-quantum, I(S:E_j) = 1 bit, and the bit is the pointer Z_S written in j's equator (the ZY channel). Even dressers are admitted: they kill the double coherences while forgiving j's record, and deg(S) ≥ 2 kills V_S. Without the direct bond there is no writer: D = ∅ with all-even P and no S–j bond is exactly dark.
- **Bell record**: Q = ∅, m ≥ 1, every dresser ratio r_k (k ∈ D) an odd integer, every private ratio r_k (k ∈ P) an even integer. Then ρ_Sj = ½|Φ^{σ₁}⟩⟨Φ^{σ₁}| + ½|Ψ^{σ₂}⟩⟨Ψ^{σ₂}| (eigenvalues ½, ½, 0, 0), I(S:E_j) = 1 bit, both marginals maximally mixed, and j's two Z_S-conditional states are identical: zero pointer content. (m ≥ 1 already kills V_S; the Bell branch does not need deg(S) ≥ 2, and K₂,₁ with deg(S) = 1 is its smallest instance.)

The exclusivity is watcher arithmetic: the pointer family needs j's shared watchers all even, the Bell family needs them all odd with m ≥ 1, and no dresser is both; for D = ∅ only the pointer family is available (given the write bond). Integer violations are exact kills: a mixed-parity D, an odd private watcher, a nonempty Q (for the Bell record), or a missing write bond (for the pointer record) each set the corresponding record to 0 identically. Non-integer ratios give the generic in-between (the forced/free pattern: parities forced, in-betweens generic), with one gated exception at the fully-shared m = 1 corner where the bit stays forced and only its channel rotates (see Edge cases).

**(ii) The letter is the dresser parity.** In the Bell case the luminous channel is decided by m mod 2:

    m odd  →  ⟨Y_S Y_j⟩ = ±1, ⟨X_S X_j⟩ = 0,
    m even →  ⟨X_S X_j⟩ = ±1, ⟨Y_S Y_j⟩ = 0,

equivalently sign(c₁c₂) = (−1)^m: the m-th shared dresser flips which Bell pair carries the bit. The direct S–j bond never moves the letter (it is invisible in both double coherences).

**(iii) The signs are explicit.** With σ₁ = sign(c₁), σ₂ = sign(c₂):

    σ₁ = Π_{k∈D} (−1)^{(1+r_k)/2} · Π_{k∈P} (−1)^{r_k/2},
    σ₂ = Π_{k∈D} (−1)^{(1−r_k)/2} · Π_{k∈P} (−1)^{r_k/2},

and the sign of the luminous correlator is σ₂ in both letter cases. On the pointer side the record's sign is Π_{k∈D∪P} (−1)^{r_k/2}, shared and private watchers alike (an r_k ≡ 2 mod 4 watcher returns the record rotated by π: ⟨Z_S Y_j⟩ = −1, whether it is a private watcher or an even dresser). This is [F135's signed-coherence corollary](PROOF_RECORD_PARITY_LAW.md) made observable: mutual information is sign-blind, the letter and the signs are where the signed cosines become measurable.

**(iv) The dephasing dressing.** With γ on, the double coherences are multiplied by e^{−2(γ_S+γ_j)t}, j's single coherences by e^{−2γ_j t}, and nothing else changes:

- the Bell record pays BOTH sites: I = 1 − h₂((1+κ)/2) with κ = e^{−2(γ_S+γ_j)t*}, and the luminous correlator reads ±κ;
- the pointer record pays only the witness (F135's Law C face): γ_S is exactly invisible to it;
- γ on dressers, watchers, or any other traced site is exactly invisible to the pair: watching the writers costs the record nothing.

## Derivation (the channel algebra of Proposition 1)

Proposition 1 factors every pair entry into the direct-bond phase, the pair's own dephasing factors, and one cosine per third site k: cos(t·[Δ_Sk·δ_S + Δ_jk·δ_j]), where δ_S = z_S − z'_S ∈ {0, ±2} and δ_j = z_j − z'_j label the channel. At t* = π/(4Δ_S), with all S-bonds equal to Δ_S, the four channel classes read:

**Populations** (δ_S = δ_j = 0): every factor 1; all four populations stay ¼ forever.

**Single-S coherences** (δ_S = ±2, δ_j = 0): every k ∈ D ∪ Q is an S-bond and contributes cos(2Δ_S t*) = cos(π/2) = 0; P contributes nothing. deg(S) ≥ 2 guarantees D ∪ Q ≠ ∅, so these vanish identically: V_S = 0, the hinge both families share.

**Single-j coherences** (δ_S = 0, δ_j = ±2): each k ∈ D ∪ P contributes cos(2Δ_jk t*) = cos(π·r_k/2); Q contributes nothing; the direct bond contributes the conditional phase e^{−it*Δ_S·z_S·δ_j} that makes this channel the pointer record. The magnitude is F135's record radius β_j = Π_{D∪P} |cos(π r_k/2)|: **perfect iff every watcher of j, shared or private, is even; zero iff any is odd.**

**Double coherences:** for c₁ (δ_S = δ_j = +2): k ∈ D contributes cos(π(1+r_k)/2), k ∈ P contributes cos(π r_k/2), k ∈ Q contributes cos(π/2) = 0. For c₂ (δ_S = +2, δ_j = −2): D contributes cos(π(1−r_k)/2), P the same cos(π r_k/2), Q again 0. The direct bond drops out of both (z_S z_j − z'_S z'_j = 0 on both channels). Hence: **Q ≠ ∅ kills both doubles; |c₁| = |c₂| = ¼ iff every dresser is odd and every private watcher even** (cos(π(1±r)/2) = ±1 exactly at odd r, 0 at even r; cos(π r/2) = ±1 at even r, 0 at odd r); the sign formulas of (iii) are these cosines read off, and sign(c₁c₂) = Π_D (−1)^{(1+r_k)/2 + (1−r_k)/2} = (−1)^m.

The mutual exclusivity in (i) is now one sentence: a perfect pointer record needs all of j's watchers even, a perfect Bell record needs the shared ones odd, and no watcher is both. A mixed D (some odd, some even) kills both families despite closure: dark.

The correlator dictionary: ⟨XX⟩ = 2(c₁+c₂), ⟨YY⟩ = 2(c₂−c₁), ⟨XY⟩ = ⟨YX⟩ ∝ Im(c) = 0 here, ⟨ZZ⟩ = 0 (populations uniform), which turns σ₁ = ±σ₂ into the letter alternation of (ii). The Bell-mixture eigenvalues {½, ½, 0, 0} and the γ-dressed I = 1 − h₂((1+κ)/2) follow by reading the two 2×2 blocks [[¼, c],[c*, ¼]].

## Corollaries

1. **The watcher parity rule.** Private watchers must be even for EITHER record to survive; the shared dressers' parity then selects the family: **all even (with the write bond) → pointer record, all odd → Bell record, mixed → dark.** The parity-by-sharedness clause (shared → odd, private → even) is the BELL family's aliveness; F135's Law A trichotomy is the D = ∅ column; the plaquette (Law B′) is the uniform m = 2 Bell row; the triangle at S, which F135's triangle-free hypothesis set aside, splits by its far bond's parity: odd → the m = 1 Bell row, even → a pointer record (both now gated).
2. **A witness never records both.** At full strength the pointer record and the Bell record exclude each other; which one a witness CAN hold is decided by the parity of its shared dressers (and the presence of the write bond), not by sharing per se.
3. **The uniform catalogue.** At uniform coupling every ratio is 1: leaves record the pointer; a triangle corner records Y⊗Y; a square's opposite corner X⊗X; K₂,m alternates YY/XX with m; every pair on the 5- and 6-cycle is exactly dark at t* (each neighbor has an odd private watcher; each non-neighbor sees a nonempty Q or an odd private watcher).
4. **The sign is an observable.** The signed-coherence corollary of F135 predicted records returning "rotated by π" under even watchers; the letter law exhibits the π: ⟨Z_S Y_j⟩ = −1 at a private r = 2 (back at +1 for r = 4), and the Bell signs walk the σ formulas (K₂,₁ at r = 3 reads YY = −1; K₂,₂ at ratios (1,3) reads XX = −1).
5. **The γ asymmetry.** The Bell record decays at the SUM of the pair's watching rates, the pointer record only at the witness's own; and no third site's watching ever touches either. In the record-price vocabulary: an anti-pointer testimony pays double, and watching the messengers is free.
6. **The fan-out: anti-pointer redundancy is not degree-bounded.** Law B bounds R_perfect by deg(S) because pointer records live on write bonds; Bell witnesses need not be neighbors. On K_{R+1,2} (S and R corners sharing the same two dressers, uniform coupling) every corner holds a perfect X⊗X bit while deg(S) = 2, and the corners Bell-record each other pairwise: a Bell-record clique through two shared dressers. The pointer's redundancy is bounded by the bonds S owns; the anti-pointer's only by who shares its dressers. (Pairwise-marginal statements; gated at K₄,₂.)

## Edge cases and scope

- **The pendant role-swap** (S's only bond is the write bond, D ∪ Q = ∅; hypothesis H excludes it, and this is why): V_S survives, the two-family classification does not apply, and the pair reads BACKWARDS. With an odd watcher on j, the pair holds exactly 1 separable bit, but it is j's pointer Z_j written in S's equator (the YZ channel): S has become the witness. The odd watcher is EXISTENTIAL: one odd-integer watcher zeroes j's radius and both double coherences at once, so further watchers, even or non-integer, do not spoil the swap (gated). Mutual information is symmetric; which qubit is "the system" is the reading's choice, and a pendant S is simply the better-protected witness of its interior neighbor. With an even watcher (or none), the pair is entangled instead: the bare dimer at θ_w = π/2 carries I = 2 bits, quantum correlation, not a classical record. Any second S-neighbor restores the hinge and closes the whole channel. (All gated; a review round caught the first hypothesis draft admitting this family unclassified.)
- **Asymmetric write couplings and J > 0** inherit [F135's scope](PROOF_RECORD_PARITY_LAW.md) unchanged: unequal S-bonds break the symmetric readout's clean endpoints, and hopping breaks the diagonal substrate entirely; everything here is the J = 0 face.
- **Non-integer ratios** are the generic in-between of (i): both letter channels partially lit, |c₁| ≠ |c₂| in general, the forced/free texture. One gated exception at the fully-shared corner: with the write bond, a SINGLE dresser and nothing else, the pair stays rank-2 and holds a full forced bit at EVERY far-bond ratio r; the ratio only rotates which channel carries it, (⟨YY⟩, ⟨ZY⟩) = (cos, sin) of π(1−r)/2, sweeping continuously between the m = 1 Bell row (r odd) and the rotated pointer record (r even). The family label varies; the bit does not.

## The frame note (what is gauge-free)

The letter names X and Y in the frame the world itself provides: the dephasing basis fixes Z, the initial state |+⟩^⊗N fixes the equatorial origin. A local Z-rotation on either site rotates equatorial operators into each other (X into Y), so the letter symbols are frame-relative, and so is the bare sign of c₁c₂. Gauge-free are: every mutual-information value, the ALTERNATION of sign(c₁c₂) with m (its ratio to the r = 1 baseline: the letter alternation itself), and all signs referenced to that baseline within one frame. The gates fix the frame once (the standard Pauli matrices in the computational basis) and read everything in it.

## The graph-state instant (where the uniform rows were already known)

At uniform coupling the readout point is special beyond our arc, and the honest statement is sharper than "a stabilizer reader will recognize this." The exact identity (checked entrywise on the bond factors):

    e^{−i(π/4)·Z_aZ_b} = CZ_ab · e^{iπ/4} · e^{−i(π/4)Z_a} · e^{−i(π/4)Z_b},

so U(t*) on |+⟩^⊗N is the graph state |G⟩ dressed by one local rotation e^{−i(π/4)·deg(l)·Z_l} per site (up to a global phase): the Ising evolution carries a DEGREE-DEPENDENT frame rotation relative to the CZ circuit.

Hein, Eisert and Briegel classified the two-qubit marginals of (unweighted) graph states in 2004 (arXiv:quant-ph/0307130; review arXiv:quant-ph/0602096): a pair marginal is maximally mixed unless the two sites' neighborhoods match, and the luminous cases are rank-2 Bell mixtures with the letter keyed to ADJACENCY: a leaf carries X⊗Z, a disconnected matched pair X⊗X, a connected matched pair Y⊗Y. Composing their classification with the degree rotation above reproduces our uniform catalogue exactly, m-parity alternation included: the Φ coherence |++⟩⟨−−| picks up the phase (π/2)·(deg_S + deg_j), which is 3π for K₂,₃ (an odd multiple of π: HEB's X⊗X turns into our Y⊗Y) and 2π for the triangle (its letter stays). This reconciliation is itself gated (undoing the degree rotation on our measured pairs recovers HEB's letters, 4 gates). So: **at uniform coupling, aliveness, the Bell-mixture form, and the letters are HEB 2004 statics seen in the physical frame; we found them independently and credit them fully.**

What the literature sweep did NOT find published, and we claim as the arc's own: (1) the integer-RATIO gating arithmetic (shared → odd, private → even, Q = ∅), the weighted generalization of "matched neighborhoods"; (2) the letter and sign closed forms in the coupling ratios, valid off the uniform point; (3) the pointer/Bell mutual exclusivity as a record statement; (4) the dephasing dressing: the both-sites price of the Bell record, the writers' exact immunity, and the two-species γ race (the rate arithmetic underneath is textbook Lindblad, the repo's own −2γk law; the two-record-species contrast is not); (5) the Darwinism packaging: the plaquette bit as a concrete, constructive realization of the Holevo/discord complementarity (Zwolak-Zurek 2013) at graph distance 2.

**Literature sweep (2026-07-18, two independent web sweeps; searches broad, not exhaustive; HEB's marginal rules verified here numerically rather than from their equations):**

- OVERLAPS (prior art for the uniform statics): Hein-Eisert-Briegel, PRA 69, 062311 (2004), arXiv:quant-ph/0307130 and the 2006 review arXiv:quant-ph/0602096 (pair marginals, Bell mixtures, adjacency letters); arXiv:1006.1594 (these marginals carry no entanglement); Zwolak-Zurek, Sci. Rep. 3, 1729 (2013), arXiv:1303.4659 (mutual information = Holevo + discord; the anti-pointer concept); spectator-site partial-trace invisibility is standard structure (also Çakmak et al., Entropy 23, 995 (2021), arXiv:2106.15629).
- ADJACENT (same objects or same flavor, not these laws): arXiv:2112.15373 (weighted-graph-state pair marginals, discord without entanglement, shared-neighbor dependence; no Bell form, no letter, no parity rule); Doucet-Deffner, PRX 14, 041064 (2024), arXiv:2405.00805 (two-body Hamiltonian classification for Darwinism at the commutation level; no ratio arithmetic); Zwolak-Riedel-Zurek, Sci. Rep. 6, 25277 (2016), arXiv:1603.01916 (spin-bath records, commensurate-time recurrences); Roszak-Korbicz, PRA 100, 062127 (2019), arXiv:1904.08261, and arXiv:1907.12447 (PRR 2020) (entanglement necessity for dephasing records); Riedel-Zurek-Zwolak, NJP 14, 083010 (2012) (redundancy timescale competition); PRB 93, 035430 (2016), arXiv:1507.04514 (even/odd bath parity controlling Ising-bath coherence plateaus); Ciampini et al., PRA 98, 020101(R) (2018) (Darwinism on photonic cluster states); arXiv:2605.15882 (complementary quantum and classical records); arXiv:2406.09956 (graph-state marginal invariants).

## Verification

All statements are gated in [`simulations/qd_letter_gates.py`](../../simulations/qd_letter_gates.py), output [`simulations/results/qd_letter/qd_letter_gates_out.txt`](../../simulations/results/qd_letter/qd_letter_gates_out.txt), **87/87 at machine zero**, every prediction fixed in the code before the run:

- letters and signs: K₂,₁ at r = 1, 3, 5 (YY = +1, −1, +1); K₂,₂ at (1,1) and (1,3) (XX = +1, −1); triangle with a private watcher at r = 2, 4 (YY = −1, +1);
- kills: even dresser, odd private watcher, nonempty Q, pentagon and hexagon (all pairs tested exactly dark);
- structure: eigenvalues {½, ½, 0, 0}, ⟨XY⟩ = ⟨YX⟩ = ⟨ZZ⟩ = 0, TD of conditionals 0 (Bell) and 1 (pointer);
- γ: I = 1 − h₂((1+κ)/2) and ⟨XX⟩ = κ at γ_S = γ_j = 0.05 (10⁻⁷); dresser-γ invisibility exact; pointer record blind to γ_S;
- the signed corollary on the pointer channel: ⟨Z_S Y_j⟩ = −1 at private r = 2, +1 at r = 4;
- exclusivity: Bell channels dead at the leaf, pointer channel dead at the plaquette;
- the pointer family with even dressers (the region a review round caught the first theorem draft mischaracterizing): the even-dresser triangle reads I = 1, ⟨Z_S Y_j⟩ = −1, TD = 1 with dead Bell channels; two even dressers read ZY = +1; the write-bond necessity (D = ∅, P even, no S–j bond → exactly dark);
- the hinge and the pendant role-swap (the region a second review round caught the hypotheses admitting unclassified): pendant S with an odd watcher reads exactly 1 separable bit in the YZ channel (j's pointer recorded by S), an even watcher reads 2 entangled bits, any second S-neighbor closes the channel;
- the fan-out: K₄,₂ pairs (0,1), (0,2), (0,3) and (1,2) all read I = 1 with ⟨XX⟩ = +1 (three Bell witnesses at deg(S) = 2, and the corners record each other: the clique);
- frame reconciliation: undoing the degree rotation on the measured pairs recovers HEB's adjacency letters (K₂,₃ X⊗X, triangle Y⊗Y, leaf X⊗Z).

The gate script never evaluates the theorem's formulas to produce its "got" values: each configuration is built as the full closed-form state (the Absorption substrate) and numerically partial-traced, so a shared error in the channel algebra cannot certify itself.

## Where this sits

Registry entry: F136 in [ANALYTICAL_FORMULAS](../ANALYTICAL_FORMULAS.md). The substrate is the [Absorption Theorem](PROOF_ABSORPTION_THEOREM.md) §4.5 through [F135's Proposition 1](PROOF_RECORD_PARITY_LAW.md); the fragment side obeys [F70's page bound](PROOF_DELTA_N_SELECTION_RULE.md). The one-four thesis's stacked-parity verdict predicted the letter's shape; the forced/free law names its texture (parities forced, in-betweens generic). The sightings and the play that found the law are recorded in the session of 2026-07-18; the arc's door catalogue lives in [experiments/QUANTUM_DARWINISM_POINTER_DOOR.md](../../experiments/QUANTUM_DARWINISM_POINTER_DOOR.md).
