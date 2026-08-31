# PROOF: Ring N=4 dihedral lock, Im_max(ring, N=4, J) = 3J, which at N=4 is (3/4)·J·N

**Status:** Tier 1 derived. The 4-cycle is the bipartite-complete graph K_{2,2}; its isotropic Heisenberg Hamiltonian factors through total sublattice spins via SU(2) Casimir, yielding four distinct levels {−2J, −J³, 0⁷, J⁵}, whose multiplicities sum to 16 = 2⁴, with max gap 3J. The Liouvillian eigenmode realising this gap is the coherence |0000⟩⟨Ψ_−| between the fully polarised state |0000⟩ (the S_z=2 member of the S_tot=2 multiplet, at E = +J) and the (S_A=1, S_B=1, S_tot=0) singlet Ψ_− at E = −2J: every basis pair (i, j) in its support has the same popcount(i⊕j) = 2, so uniform dephasing acts on the whole operator as the scalar −4γ, and λ = −4γ − 3iJ, with the eigenoperator residual at machine precision across γ = 0.005 to 50 (`simulations/ring_n4_lock_gate.py` §3). The reversed ordering carries the conjugate.
**Date:** 2026-05-19
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Distinct from:** [`PROOF_RING_GAP_DOMINANCE.md`](PROOF_RING_GAP_DOMINANCE.md), the *XY* ring gap-dominance result (max|Im| = 2J = J·ρ, the dihedral lock). Same words "ring N=4 dihedral lock", different Hamiltonian and different result; this proof is the *isotropic-Heisenberg* ring (max|Im| = 3J via the K_{2,2} Casimir gap).

## Abstract

The 4-cycle is the one ring size that coincides with the bipartite-complete graph K_{2,2}: two sublattices of two qubits, every inter-sublattice bond present, none within. On it the isotropic Heisenberg Hamiltonian factors through total sublattice spins, H = J·S_A·S_B, and the SU(2) Casimir spectrum {−2J, −J³, 0⁷, J⁵} has maximal gap 3J. The six rows of the (S_A, S_B) table in Section 3 carry only four distinct values: three of those rows sit at 0, with dimensions 3 + 1 + 3, so the zero level is seven-fold. Under uniform Z-dephasing the Liouvillian's largest imaginary eigenvalue is pinned exactly to that gap,

    Im_max(ring, N=4, J) = (3/4)·J·N = 3J,   equivalently   Im/σ = 3Q/4,

independently of γ and Q = J/γ. Because the jump operators Z_l are Hermitian, the dissipator is self-adjoint in the Hilbert-Schmidt inner product, and no eigenmode can exceed the Hamiltonian's spread. That the bound is reached rather than approached takes two further hypotheses, uniform γ across the full single-site jump set and a fully polarised H extreme, both discharged in Section 4 and fenced in "Scope".

What is the 4-cycle's own here is the closed form 3J, and nothing else in the sentence. The equality max|Im λ_L| = ΔE_max(H) is [F148](../ANALYTICAL_FORMULAS.md#f148-the-imaginary-reach-is-the-hamiltonian-spread-on-every-graph-minted-2026-07-31), which holds on every graph, so the Liouvillian half of this proof is not N=4's and not the ring's. Nor does the *value* 3J single the ring out: it is the maximum of ΔE_max over connected graphs at N=4, and ten of the 38 connected labelled graphs attain it, among them the complete K_4, which is not bipartite and whose Casimir multiplet structure {−1.5², −0.5⁹, +1.5⁵} is entirely different. Bipartite-completeness therefore cannot be the cause of the number; it is one route to computing it. (The historical name "dihedral lock" is a stable identifier inherited from the XY sibling; the operative mechanism here is the SU(2) Casimir multiplet structure on K_{2,2}, not the dihedral point group, D₄'s irreps cap at dimension 2 and cannot produce the 3-fold/5-fold degeneracies that set the 3J gap.) What does not carry to larger even rings is the *rational* closed form: they keep a Q-universal lock but at irrational algebraic constants (ring N=6 at 0.7171·J·N = ((5+√13)/12)·J·N, descending toward ln 2), because the bipartite-complete structure is special to N=4. Typed as RingN4DihedralLockClaim.

## Statement

For the open quantum system on N=4 qubits with

- Hamiltonian: isotropic Heisenberg H = (J/4) Σ_{(i,j)∈E} (X_i X_j + Y_i Y_j + Z_i Z_j) on the 4-cycle bonds E = {(0,1), (1,2), (2,3), (3,0)};
- Dissipation: uniform Z-dephasing γ per site;

the Liouvillian L = −i[H, ·] + D[Z_l] satisfies the saturation

    Im_max(ring, N=4, J)  ≡  max_{λ ∈ σ(L)} |Im(λ)|  =  (3/4) · J · N  =  3 · J         (for N=4)

independently of γ and of the corresponding dimensionless ratio Q = J/γ. The equivalent dimensionless statement is

    Im_max / σ  =  3Q/4              with σ = N·γ.

## Empirical anchors (6 Q-values × γ₀=0.05, machine precision)

Q-sweep on 2026-05-19 (`simulations/f1_q_sweep_anchor.py`, output under `simulations/results/q_sweep_anchor/ring_N4_Q*.json`):

| Q | predicted 3Q/4 | observed Im/σ | rel. error |
|---:|---:|---:|---:|
| 0.5    | 0.375000 | 0.375000 | 1.2e-15 |
| 1.0    | 0.750000 | 0.750000 | 1.5e-16 |
| 1.5    | 1.125000 | 1.125000 | 7.9e-16 |
| √3     | 1.299038 | 1.299038 | 1.0e-15 |
| 2.0    | 1.500000 | 1.500000 | 3.3e-15 |
| 2.5    | 1.875000 | 1.875000 | 5.1e-15 |

All six anchors hit the prediction to within machine precision (relative error < 1e-14). The lock is Q-universal: the absolute Im_max value scales with J, but the dimensionless `Im/σ = 3Q/4` ratio is exact at every Q.

The observed column is the two-sided `max|Im λ|` the Statement is about, read off the stored sweep as `max(MaxImag, |MinImag|)`; at Q = 0.5, 1.0 and √3 it is `|MinImag|` that is marginally the larger, which is why those three error entries differ from the sweep's one-sided `MaxImag` reading. The distinction lives in the last bit and moves nothing else. `simulations/ring_n4_lock_gate.py` §8 rebuilds both columns from the JSON and checks them against the values the typed claim stores.

## Proof

### Section 1. The 4-cycle is the bipartite-complete graph K_{2,2}

The ring on N=4 sites has bond set E = {(0,1), (1,2), (2,3), (3,0)} (the 4-cycle). Partition the sites into sublattices A = {0, 2} and B = {1, 3} (two-colouring of the bipartite 4-cycle). The bond set then reads E = {(0,1), (0,3), (2,1), (2,3)}: every pair (a, b) with a ∈ A and b ∈ B is a bond, and no within-sublattice bonds exist. This is exactly the bipartite-complete graph K_{2,2} on 2+2 sites.

The 4-cycle C_4 and the bipartite-complete K_{2,2} are isomorphic as graphs. This is a coincidence specific to N=4: for N=6 the 6-cycle C_6 is bipartite (sublattices {0,2,4} and {1,3,5}) but has only 6 bonds versus K_{3,3}'s 9, so the bipartite-complete decomposition does not apply.

### Section 2. Total-sublattice-spin factorisation

Let S_i = (S^x_i, S^y_i, S^z_i) be the spin-1/2 operator on site i. Vector operators carry no arrow here: S_i and the sublattice totals S_A, S_B are three-component operators and the dot is their scalar product, while the same letters read as quantum numbers wherever they carry no dot and no square, as in S_A ∈ {0, 1} or S_tot(S_tot+1) in Section 3. Using S_i · S_j = (1/4)(X_i X_j + Y_i Y_j + Z_i Z_j), the bond Hamiltonian is J · S_i · S_j and the total Hamiltonian on K_{2,2} is

    H  =  J · Σ_{a ∈ A, b ∈ B} S_a · S_b
       =  J · (S_0 + S_2) · (S_1 + S_3)
       =  J · S_A · S_B

with total sublattice spins S_A := S_0 + S_2 and S_B := S_1 + S_3. The Hamiltonian is bilinear in the sublattice total spins, with no within-sublattice term.

### Section 3. Casimir spectrum

Using the standard total-spin Casimir identity S_A · S_B = (1/2)(S²_tot − S²_A − S²_B) with S_tot := S_A + S_B,

    H  =  (J/2) · (S²_tot − S²_A − S²_B).

Each sublattice contains two spin-1/2 sites, so S_A, S_B ∈ {0, 1}. Each (S_A, S_B) pair has its own internal Clebsch-Gordan multiplicity m_inner from coupling 1/2 ⊗ 1/2 → S_A (and 1/2 ⊗ 1/2 → S_B): m_inner = 1 for the singlet S = 0, m_inner = 1 for the triplet S = 1 (there is exactly one way to make each from two spin-1/2). Then coupling S_A and S_B gives S_tot ∈ {|S_A − S_B|, ..., S_A + S_B}, with the (2S_tot + 1) M_tot states inside each S_tot multiplet. The full eigenvalue list and the dimensions are:

| (S_A, S_B) | S_tot | E = (J/2)·(S_tot(S_tot+1) − S_A(S_A+1) − S_B(S_B+1)) | dimension = m_inner(A) · m_inner(B) · (2S_tot+1) |
|---|---:|---:|---:|
| (0, 0) | 0 | (J/2)·(0 − 0 − 0) = **0**     | 1 · 1 · 1 = 1 |
| (0, 1) | 1 | (J/2)·(2 − 0 − 2) = **0**     | 1 · 1 · 3 = 3 |
| (1, 0) | 1 | (J/2)·(2 − 2 − 0) = **0**     | 1 · 1 · 3 = 3 |
| (1, 1) | 0 | (J/2)·(0 − 2 − 2) = **−2J**   | 1 · 1 · 1 = 1 (outer-CG singlet of two inner triplets) |
| (1, 1) | 1 | (J/2)·(2 − 2 − 2) = **−J**    | 1 · 1 · 3 = 3 (outer-CG triplet) |
| (1, 1) | 2 | (J/2)·(6 − 2 − 2) = **+J**    | 1 · 1 · 5 = 5 (outer-CG quintuplet) |

The inner-CG multiplicities m_inner(A) and m_inner(B) are both 1 because two spin-1/2's couple into each of S = 0, S = 1 in exactly one way. The dimensions tracked are then just (2S_tot + 1) per row. Total Hilbert space dimension: 1 + 3 + 3 + 1 + 3 + 5 = 16 = 2^4 ✓.

The eigenvalue multiset of H is therefore

    σ(H)  =  { −2J, −J·(triplet, mult 3), 0·(seven-fold), +J·(quintuplet, mult 5) }.

Maximum eigenvalue is +J (the ferromagnetic S_tot = 2 multiplet); minimum is −2J (the perfect singlet of two anti-aligned triplet dimers). Maximum H eigenvalue gap is therefore

    ΔE_max(H_K22)  =  E_max − E_min  =  J − (−2J)  =  3J  =  (3/4) · J · N         for N=4.

### Section 4. Liouvillian eigenmode realising the bound

The Lindblad Liouvillian L = −i[H, ·] + D where D[ρ] = Σ_l γ_l (Z_l ρ Z_l − ρ) is the pure-dephasing dissipator. For any pair of H-eigenstates |α⟩, |β⟩ with eigenvalues ω_α, ω_β, the rank-1 operator |α⟩⟨β| is an eigenoperator of `−i[H, ·]` with eigenvalue `−i(ω_α − ω_β)` (so Im(λ_L) = −(ω_α − ω_β)). D is self-adjoint and negative-semidefinite in the operator inner product, which is what gives the bound of Section 5.

D acting as a **scalar** on such a rank-1 operator is a further condition, and it is not automatic: D is diagonal in the computational coherence basis, not in the H-eigenbasis, so for general |α⟩, |β⟩ it mixes these products and moves the frequency. Since D acts on |i⟩⟨j| by the rate −2γ·popcount(i⊕j), the condition is that **popcount(i⊕j) be constant across the support** of |α⟩⟨β| (see [PROOF_STAR_OPTICAL_CONFOCAL_SATURATION](PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md) Section 6 for the general form of this step).

Definite Hamming weight of each factor separately is not enough, and at N=4 it fails: |Ψ_−⟩ has definite weight 2, and so does the M=0 member of the quintuplet, yet their coherence spreads over popcount(i⊕j) ∈ {0, 2, 4} and carries three different rates at once.

The operator that does satisfy the condition is

    |Ψ_+⟩  =  |0000⟩,  the fully polarised S_z = +2 member of the S_tot = 2 multiplet (E = +J),
    |Ψ_−⟩  =  the (S_A=1, S_B=1, S_tot=0) singlet (E = −2J), of definite Hamming weight 2,

because a ket that is a single computational state forces popcount(i⊕j) = 2 throughout. Then |0000⟩⟨Ψ_−| is a Liouvillian eigenoperator with

    λ  =  −2γ·2 − i · (J − (−2J))  =  −4γ − 3i·J,        |Im λ| = 3J,

matching the empirical Im_max = (3/4)·J·N. The same holds for |1111⟩, and the reversed orderings give the conjugates. **The other three members of the quintuplet do not work**: measured eigenoperator residuals at γ = 0.5 are 1.000 (M = ±1) and 1.155 = 2/√3 (M = 0) against machine precision for M = ±2, so the failure is structural and not numerical. Note that all five have the same Rayleigh quotient ⟨v, Lv⟩ = −4γ − 3iJ, so a Rayleigh-quotient check does not see the difference; the residual does.

### Section 5. No mode exceeds the bound

Write L = A + D on operator space with the Hilbert-Schmidt inner product, A = −i[H, ·] and D the dephasing dissipator. A is skew-adjoint with spectrum {−i(ω_α − ω_β)}, so ‖A‖ = ΔE_max(H); D is self-adjoint. For a unit right eigenvector v of L with eigenvalue λ, λ = ⟨v, Lv⟩ = ⟨v, Av⟩ + ⟨v, Dv⟩, and the second term is real. Hence Im(λ) = Im⟨v, Av⟩, and the maximum |Im(λ_L)| over the L-spectrum is bounded by

    max |Im(λ_L)|  ≤  max{|ω_α − ω_β| : ω_α, ω_β ∈ σ(H)}  =  ΔE_max(H).

For K_{2,2} this is 3J. Combined with the realising mode in Section 4, the bound is achieved exactly:

    Im_max(ring, N=4, J)  =  ΔE_max(H_K22)  =  3J  =  (3/4) · J · N.

The same bound argument applies under non-uniform γ_l per site, as long as the dissipator is pure-dephasing (Z_l jump operators only): only self-adjointness of D was used, not the value of the rates. Read the bound for what it is, an interval and not a fixed point: individual frequencies do move under the watching, and most of them leave the difference set {ω_α − ω_β} (at N = 4, J = 1, γ = 0.5, 28 of the 35 distinct signed imaginary parts are not in it, counting values as distinct at a separation of 1e-6). What no eigenmode can do is leave the interval.

**Only the bound survives non-uniform rates, not the saturation.** Section 4's realiser needs D to act on the coherence as one scalar, and site-dependent γ_l destroys that. Measured at J = 1: γ = (1, 0.5, 0.5, 0.5) gives 2.7247448714 and γ = (1, 0.1, 1, 0.1) gives 2.3285347127, both strictly below 3J and both under the bound. Setting a rate to zero is the same break seen from the other side, since it removes that site's jump operator altogether: γ = (2, 0, 0, 0) is Z_0 alone at rate 2, and gives 2.3827292861. Even the full jump set is not enough if it is not the single-site one: the two-body Z_0Z_1 at rate 0.5 gives 2.9433995088, and Z_0 alone at that rate gives 2.8961214544. All of these are rate-dependent, as they must be once the mode is no longer an eigenoperator: Z_0 alone climbs to 2.9990042309 at rate 0.05. The hypothesis the equality needs is the **full single-site set {Z_l} at one common rate**, and only there is the value γ-independent.

### Section 6. Q-universality

The formula Im_max = (3/4)·J·N depends on J but not on γ. Translating into the dimensionless ratio Im/σ where σ = N·γ:

    Im_max / σ  =  (3/4) · J · N / (N · γ)  =  (3/4) · (J/γ)  =  (3/4) · Q.

This is the Q-universal lock observed in the Q-sweep table (Section "Empirical anchors").

## Scope: four hypotheses, each with a fence

Everything above is stated at uniform γ with the full single-site jump set. Each of the four hypotheses below is load-bearing, and each fence is a measurement in `simulations/ring_n4_lock_gate.py`.

- **Uniform γ across all N sites.** Required for the equality; the bound survives without it. γ = (1, 0.5, 0.5, 0.5) gives 2.7247448714 against 3J (§5).
- **The full single-site set {Z_l}.** Required for the equality. At a common rate 0.5, Z_0 alone gives 2.8961214544 and the two-body Z_0Z_1 gives 2.9433995088, both Hermitian. Neither value is γ-universal, which is the same break read off the other axis.
- **A fully polarised H extreme.** Required. This, not isotropy, is what carries the attainment: XXZ at Δ = 2 saturates its own spread 4.7320508076 = 3 + √3 exactly at γ = 0.05, 0.5 and 2, Δ = −1 saturates at 3.0000000000, pure Ising at 2.0000000000, while Δ = 0.5 fails (2.8672194968 at γ=0.05 down to 2.1861406616 at γ=2, against a spread of 2.8722813233) and the transverse-field Ising fails. The dividing line is exactly whether |0000⟩ or |1111⟩ is an extremal eigenvector of H, which for XXZ means |Δ| ≥ 1. Saturation and γ-independence stand or fall together: every non-saturating row is also γ-dependent.
- **Hermitian jump operators.** Required, or the **bound itself** fails. H = 0 on one qubit with c = I + iY has spectrum {0, 0, −2 ± 2i}, so max|Im| = 2 against a bound of 0. The informal justification that the dissipator "adds only real decay" would have licensed the bound there too; the load-bearing fact is that D is Hilbert-Schmidt self-adjoint, and there ‖D − D†‖_F = 4√2 = 5.6569.

## Why this is N=4-specific

The bipartite-complete structure C_4 = K_{2,2} relies on the 4-cycle having exactly 4 bonds (one per (A, B) pair). For longer cycles:

- 6-cycle has 6 bonds, K_{3,3} has 9: a 6-cycle is bipartite but NOT bipartite-complete.
- 8-cycle has 8 bonds, K_{4,4} has 16: same story.

For odd N (3-cycle, 5-cycle, ...) the ring is not even bipartite. So the K_{2,2} = C_4 coincidence is unique to N=4. Ring N=6, ring N=8 etc. show Q-universal locks too (empirically, [`hypotheses/F1_DISSIPATION_GAP_PATTERN.md`](../../hypotheses/F1_DISSIPATION_GAP_PATTERN.md): ring N=6 = 0.717129·J·N at 6 Q-anchors), but the per-N constant is no longer rational, and no Casimir argument as simple as the N=4 one produces it. It does stay elementary for a while: 4H has integer entries, so E₀(N) is an algebraic number, and factoring the exact characteristic polynomial of the S_z = 0 sector over ℚ gives c₆ = (5+√13)/12 (the ground-state factor is λ² + 8λ − 36 **in units of 4H**, i.e. 12c² − 10c + 1 in c itself), c₈ the largest root of 512c³ − 640c² + 232c − 25, and c₁₀ the largest root of a sextic (`ring_dihedral_lock_limit.py` STAGE 2 checks each c_N against its minimal polynomial in c; the 4H-scaled factors are not printed there). What needs Bethe ansatz is the N → ∞ limit, not the individual N.

## The N → ∞ limit: c_∞ = ln 2 (resolved 2026-06-04)

Although the per-N value needs Bethe ansatz, the LIMIT is closed. Section 5 gives the bound Im_max(L) ≤ ΔE_max(H) at every N (the jump operators are Hermitian, so D is self-adjoint in the operator inner product), and what turns the bound into the equality Im_max(L) = ΔE_max(H) = E_max − E_min, reducing the lock to the Hamiltonian alone, is [F148](../ANALYTICAL_FORMULAS.md#f148-the-imaginary-reach-is-the-hamiltonian-spread-on-every-graph-minted-2026-07-31), at every N and on every graph. Section 4's realiser is the N=4 instance of F148's general one: what it needs is that |0…0⟩ be a Z-product state and an H extreme, and that [H, Σ_l Z_l] = 0 put an E_min eigenvector inside one Hamming rung k*, so that uniform dephasing acts on |0…0⟩⟨β| as the scalar −2γk*. Nothing in that is the 4-cycle's or N=4's. So the equality holds at every N of the table below, odd rows included, with no assumption carried and no 4^N diagonalisation needed: `ring_dihedral_lock_limit.py` STAGE 1(a) certifies it at N = 4…16 even and 5…13 odd, `ring_n4_lock_gate.py` §7b at N = 6, 8, 10. On that footing the dimensionless constant is

    c_N  ≡  Im_max / (J·N)  =  (E_max − E_min) / (J·N)  =  1/4 − E₀(N)/(J·N),

with E_max = J·N/4 (the ferromagnet, exact) and E₀(N) the antiferromagnetic Heisenberg-ring ground state. The per-bond ground energy of the spin-½ Heisenberg ring has the Bethe/Hulthén thermodynamic limit E₀/(J·N) → 1/4 − ln 2, hence

    c_∞  =  1/4 − (1/4 − ln 2)  =  ln 2  =  0.693147…     (NOT 1/√2 = 0.707107).

Computing E₀(N) directly ([`simulations/ring_dihedral_lock_limit.py`](../../simulations/ring_dihedral_lock_limit.py), rung-resolved ground state, N = 4..16) confirms it:

| N | c_N = 1/4 − E₀/N | c_N − ln 2 |
|---:|---:|---:|
| 4 | 0.75000 | +0.05685 |
| 6 | 0.71713 | +0.02398 |
| 8 | 0.70639 | +0.01324 |
| 10 | 0.70154 | +0.00840 |
| 12 | 0.69895 | +0.00580 |
| 14 | 0.69740 | +0.00425 |
| 16 | 0.69639 | +0.00325 |

Two independent Liouvillian-side measurements meet the Hamiltonian-side table where they overlap, which is a check on the reduction rather than a premise for it. The N=6 value 0.71713 reproduces the F1_DISSIPATION_GAP empirical 0.717129. At N=8 the agreement is sharper than any document has been carrying: the F1 SLOW_N8 sweep of 2026-05-18 formed the whole 4⁸ spectrum once and recorded `MaxImag = 5.651093408937174` (`simulations/results/f1_n8_n9_metrics/ring_N8.json`, J=1, γ=0.5), against the EXACT value, which needs no second eigensolver: ΔE_max(8) = 8·c₈, and c = x/8 clears the cubic below to the monic integer `Q(x) = x³ − 10x² + 29x − 25`, whose largest root it is. Evaluated over the rationals, Q places the archived double exactly, two doubles below the true algebraic number and nowhere else, while the published form of that same datum has been the four-digit `Im/σ = 1.4128 = 0.7064·Q`. Along the even N the sequence decreases to ln 2 at the rate the periodic c = 1, v = π/2 conformal finite-size form predicts, E₀(N)/N = e_∞ − π²/(12N²), so N²(c_N − ln 2) falls monotonically toward π²/12 = 0.822467: measured 0.909645, 0.863355, 0.847328, 0.839745, 0.835497, 0.832846, 0.831064 at N = 4…16, the remaining 1.0% at N=16 being the chain's known logarithmic correction. That law is what `ring_dihedral_lock_limit.py` STAGE 0 gates, and the refuted reading fails it outright: N²(c_N − 1/√2) changes sign at N=8 and diverges. It passes 1/√2 between N=6 and N=8 (c₆ = 0.717129 above, c₈ = 0.706387 below, and no N sits on it); 1/√2 is only a value the sequence steps across, not the limit (the same red-herring lesson as the birth-canal s* = 0.709). So the Q-universal ring dihedral lock, left open for general N, has the exact limit c_∞ = ln 2.

**The table is the even branch, and the sequence over all N is not monotone.** The odd rings are frustrated, their antiferromagnetic ground state sits higher per site, and their per-site spread is correspondingly smaller: c₅ = 0.623607, c₇ = 0.657883, c₉ = 0.671922, c₁₁ = 0.678994. Those lie **below** ln 2 and rise toward it while the even ones fall toward it. Two monotone branches closing on the same limit from opposite sides, not one descent:

| N | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|
| c_N | 0.750000 | 0.623607 | 0.717129 | 0.657883 | 0.706387 | 0.671922 | 0.701545 | 0.678994 | 0.698949 |

The even rows alone do not determine the shape of the sequence. The neighbouring parity does, and it runs the other way.

## Relationship to the star saturation

[STAR_CONFOCAL_LIMIT.md](../../experiments/STAR_CONFOCAL_LIMIT.md) shows the analogous result for the star topology:

    Im_max(star, N, J)  =  J·N/2

via a parallel hub-spoke Casimir construction H = J · S_0 · S_L with S_L = Σ leaf spins. The Casimir spectrum has max gap J·N/2 (the maximally-ferromagnetic-leaves S_L = (N-1)/2 sector flipping the hub). The ring N=4 result is structurally the same kind of object: a topology where the SU(2)-invariant Heisenberg Hamiltonian factors through a sum-of-Casimirs form, so that ΔE_max has a closed form. What turns ΔE_max into Im_max is F148 in both cases, not the Casimir.

The ratio between the two is `(3/4)·J·N / (J·N/2) = 3/2`: ring N=4 carries 50% more imaginary spread than star N=4. Both sit at an end of the same range. J·N/2 is the **minimum** of ΔE_max over connected graphs, with the star its unique minimiser ([F147](../ANALYTICAL_FORMULAS.md#f147-the-star-spread-the-hub-leaf-casimir-gives-jn2-the-smallest-hamiltonian-spread-among-connected-graphs-at-n--6-derived-2026-05-19-registered-2026-07-31), exhaustive at N ≤ 6); at N=4 the 4-cycle sits at the **maximum** of the same quantity. The two ends are 2J and 3J, and the ratio is that of the two ends, not a property either graph carries alone.

The maximum, unlike the minimum, has no unique attainer and does not survive the next N. Ten of the 38 connected labelled graphs on four vertices reach 3J: the three labellings of C_4, the six of the diamond K_4 − e, and K_4. And at N=5 the ring is not even close to the top, 3.1180339887 against K_5's 4.0000000000; at N=6, 4.3027756377 against 6.0000000000. So the coincidence "the 4-cycle is a maximiser" belongs to N=4 alone, exactly like the rational value.

## Verification

- Gate: [`simulations/ring_n4_lock_gate.py`](../../simulations/ring_n4_lock_gate.py). Every N=4 number this document states is produced there: the Casimir table with its multiplicities, the lock over J and over four decades of γ, the six Q-anchors against the stored JSON, the realiser and the three quintuplet members that fail it, the 10-of-38 tie at 3J, the four scope fences, and the 28-of-35 count. It carries the ladder's low rungs too, c₄/c₆/c₈ even and c₅/c₇/c₉/c₁₁ odd, and the F148 certificate at N = 6, 8, 10. The rows from N = 10 up, c₁₀'s sextic, the π²/12 law and the N=8 archive comparison belong to the sibling named in the next bullet.
- Python anchor at 6 Q-values × γ₀=0.05: [`simulations/f1_q_sweep_anchor.py`](../../simulations/f1_q_sweep_anchor.py) → `simulations/results/q_sweep_anchor/ring_N4_Q*.json`.
- Typed claim: [`compute/RCPsiSquared.Core/Symmetry/RingN4DihedralLockClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/RingN4DihedralLockClaim.cs) (Tier 1 derived) with `Predict(J)` returning `(3/4) · J · N = 3J` at N=4.
- The larger rings: [`simulations/ring_dihedral_lock_limit.py`](../../simulations/ring_dihedral_lock_limit.py). STAGE 0 the c_N ladder, both parity branches gated, and the ln 2 limit gated as the π²/12 finite-size law rather than as a threshold; STAGE 1 F148 from below, (a) the |ferro⟩⟨ground| attainment certificate at every N of the ladder with no 4^N object formed, (b) the bound cross-checked on the full Liouvillian spectrum at N=4 and, under `--slow`, N=6, both at two γ, (c) the archived N=8 full-spectrum datum; STAGE 2 the finite-N closed forms against their integer minimal polynomials, N = 5, 6, 8 and 10, the odd rung among them showing the odd branch is algebraic too.

## Cross-references

- Parent: [F1PalindromeIdentity](../../compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs) (the F1 master under which this Im-max bound is verified by the same SLOW_N* sweep infrastructure).
- Sister Im-max bound, closed in the same style on 2026-05-19: [the star optical-confocal saturation proof](PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md), with the record in [STAR_CONFOCAL_LIMIT.md](../../experiments/STAR_CONFOCAL_LIMIT.md).
- Sister Q-universal lock: ring N=6 at 0.717129·J·N = ((5+√13)/12)·J·N, a quadratic surd from the exact S_z=0 characteristic polynomial (see [F1_DISSIPATION_GAP_PATTERN.md](../../hypotheses/F1_DISSIPATION_GAP_PATTERN.md) "Ring N=6 dihedral lock" section). The constant that stays open is the general-N one.
- Companion typed claim from the same May 2026 sharpening sprint: [F4KernelDimensionByComponentsClaim](../../compute/RCPsiSquared.Core/Symmetry/F4KernelDimensionByComponentsClaim.cs) (kernel-dim factorisation across components, Tier 1 derived 2026-05-19).
- Q-anchor canonical table: [`docs/Q_REGIME_ANCHORS.md`](../Q_REGIME_ANCHORS.md).
