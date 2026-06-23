# PROOF: F49 cross-term `‖{L_H, L_Dc}‖²` closed form under non-uniform γ

**Status:** Tier 1 derived. Extends the F49 cross-term closed form `‖{L_H, L_Dc}‖² = 4γ²·(N−2)·‖L_H‖²` (uniform γ) to arbitrary site-dependent {γ_l}. The non-uniform formula splits into a spectator part (per-bond, depending on Σ_{m∉bond} γ_m²) plus a bond-asymmetry part (per-bond, depending on (γ_i − γ_j)² with a per-Pauli-class coefficient G(bond, H)). The bond-asymmetry coefficient G(bond, H) is exactly 4·‖L_{ZZ-class part of H}^bond‖²; for canonical H-classes this gives (4/3)·‖L_H^bond‖² (Heisenberg), 4·‖L_H^bond‖² (Ising), 0 (XY and soft XY+YX). Bit-exact verification at N = 3, 4, 5 across all four H classes.
**Date:** 2026-05-18
**Authors:** Thomas Wicht, Claude (Opus 4.7)

## Abstract

F49 closed the cross-term Frobenius norm `‖{L_H, L_Dc}‖²` for uniform Z-dephasing: a clean per-bond product of the Hamiltonian commutator norm and a γ² factor times (N − 2). This proof extends it to arbitrary site-dependent γ rates. The natural extension would replace the (N − 2) factor by something involving Σ γ², but the truth is more interesting.

The non-uniform closed form splits into two structurally distinct pieces. The first is a spectator part: each bond contributes a term proportional to its bond-local commutator norm times the sum of squared γ rates AT THE NON-BOND SITES (the sites the bond does not touch). The second is a bond-asymmetry part: each bond contributes a term proportional to the SQUARED DIFFERENCE between its two endpoint γ rates, times a per-bond coefficient that depends on the Pauli-class composition of the bond's Hamiltonian.

The bond-asymmetry coefficient G(bond, H) is the diagnostic gem. It works out to exactly four times the Frobenius norm of the "ZZ-class part" of the bond Hamiltonian. For canonical bond classes the coefficient is clean: Heisenberg picks up (4/3) of the bond norm, Ising picks up the full 4 (since ZZ is the entire bond), XY picks up zero (no ZZ content), and the F87-soft XY+YX combination also picks up zero. The non-uniform γ generalization is therefore not just a numerical extension; it gives a per-bond fingerprint of how much of each bond's Hamiltonian sits in the ZZ direction.

The diagnostic upshot is that measuring the F49 cross-term norm under non-uniform γ reads off both the bond Hamiltonian composition (via the asymmetry coefficient G) and the spectator γ rates (via the spectator part). For uniform γ the asymmetry term vanishes and the F49 closed form is recovered. The non-uniform extension is therefore a strict refinement, not a replacement, and the ZZ-content per bond becomes measurable from outside the system through the cross-term observable.

## Statement

Let H be a 2-bilinear bond Hamiltonian on N qubits, `H = Σ_b H_b` with each bond term `H_b = Σ_t c_t · α_t^{(i_b)} β_t^{(j_b)}` for Paulis α_t, β_t ∈ {X, Y, Z} on bond sites (i_b, j_b). Let L_D be the site-dependent Z-dephasing dissipator with rates {γ_l},

    L_D(ρ) = Σ_l γ_l · (Z_l ρ Z_l − ρ),    {γ_l ≥ 0}_l=1^N arbitrary,

let σ := Σ_l γ_l, and define the F1-centered dissipator L_Dc := L_D + σ · I.

**Theorem (non-uniform γ cross-term closed form).** The F1 cross-term Frobenius norm satisfies

    ‖{L_H, L_Dc}‖²_F
       =  4 · Σ_b ‖L_H^bond_b‖²_F · Σ_{m ∉ bond_b} γ_m²        (spectator part)
        + Σ_b G(bond_b, H) · (γ_{i_b} − γ_{j_b})²                (bond-asymmetry part)

where `‖L_H^bond_b‖²_F` is the Frobenius norm squared of the per-bond Hamiltonian superoperator (built on the full N-qubit Pauli-string operator space; for a single Heisenberg bond at J = 1 this equals 384 at N = 3, scaling as 4^(N−2) · 96 in general; see [Conventions](#conventions)), and the per-bond per-class bond-asymmetry coefficient is

    G(bond_b, H)  =  4 · ‖L_{ZZ-class part of H_b}^bond_b‖²_F.

For canonical H-classes, the ZZ-fraction of `‖L_H^bond_b‖²_F` is determined by the bond Pauli decomposition, giving:

| H class                | G(bond_b, H) / ‖L_H^bond_b‖²_F | A-classification reading |
|------------------------|--------------------------------|--------------------------|
| Heisenberg J·(XX+YY+ZZ) | **4 / 3**                     | ZZ contributes 1/3 of the bond norm; XX, YY contribute 2/3 with A = 0 |
| Ising J·ZZ              | **4**                         | ZZ contributes 100% of the bond norm |
| XY J·(XX+YY)            | **0**                         | no ZZ component; XX, YY both have A = 0 |
| Soft Π²-odd J·(XY+YX)   | **0**                         | no ZZ component; XY, YX both have A = 0 |

The uniform-γ identity ([F49 / PROOF_CROSS_TERM_FORMULA](PROOF_CROSS_TERM_FORMULA.md))

    ‖{L_H, L_Dc}‖²_F  =  4γ² · (N−2) · ‖L_H‖²_F

is recovered as the special case γ_l ≡ γ (the bond-asymmetry part vanishes because γ_i = γ_j on every bond; the spectator part collapses to `4γ² · (N−2) · Σ_b ‖L_H^bond_b‖²_F = 4γ² · (N−2) · ‖L_H‖²_F` by [PROOF_CROSS_TERM_FORMULA Step 4 disjoint-bond-supports lemma](PROOF_CROSS_TERM_FORMULA.md)).

## Conventions

- **Pauli letters** (a, b) ∈ {(0,0), (1,0), (0,1), (1,1)} = (I, X, Z, Y) following the framework's Klein-Vierergruppe convention ([`framework/pauli.py`](../../simulations/framework/pauli.py)). `bit_a(α) = 1` iff α ∈ {X, Y}; `bit_b(α) = 1` iff α ∈ {Y, Z}.
- **Pauli-string basis on N sites** is the 4^N orthonormal basis {σ_α} with Tr(σ_α^† σ_β) / 2^N = δ_{αβ}. The Frobenius norm of a 4^N × 4^N superoperator is computed in this orthonormal basis (Hilbert-Schmidt product).
- **Hamiltonian superoperator** L_H = −i [H, ·] in the Pauli-string basis; per-bond L_H^bond_b = −i [H_b, ·] built on the full N-qubit Pauli-string operator space, including spectator I-tensors. Concretely, `‖L_H^bond_b‖²_F` for a single Heisenberg J = 1 bond equals 384 at N = 3, 1536 at N = 4, 6144 at N = 5 (each step multiplies by 4 because each additional spectator site contributes a factor `tr(I_4) = 4` in the Frobenius-norm tensor calculation; intrinsic local 2-qubit Heisenberg-bond superoperator norm is 96). The closed form here uses this convention throughout, matching `_bond_LH_norm_sq(N, bond, terms)` in [`simulations/f49_nonuniform_gamma_crossterm_verify.py`](../../simulations/f49_nonuniform_gamma_crossterm_verify.py).
- **Z-dephasing dissipator** in the Pauli basis acts as a diagonal multiplier on each σ_α with eigenvalue `d_α = −2 Σ_l γ_l · bit_a(α_l)`, the per-site weighted variant of the F49 uniform `d_α = γ · (N − 2 · w_XY(α))` ([PROOF_CROSS_TERM_FORMULA Step 1](PROOF_CROSS_TERM_FORMULA.md)). The F1-centered dissipator L_Dc := L_D + σ · I is also diagonal, with eigenvalue `d_α^c := σ − 2 Σ_l γ_l · bit_a(α_l) = Σ_l γ_l · (1 − 2·bit_a(α_l)) = Σ_l γ_l · ε_l(α)` where ε_l(α) := 1 − 2·bit_a(α_l) ∈ {+1, −1}. (Equivalently, ε_l(α) = +1 if α_l ∈ {I, Z} and −1 if α_l ∈ {X, Y}; matches the F49-proof ε notation.)

## Empirical anchor

Phase 1 verification at commit `1c6701c` ([`simulations/f49_nonuniform_gamma_crossterm_verify.py`](../../simulations/f49_nonuniform_gamma_crossterm_verify.py)) confirmed the formula bit-exact at N = 3, 4, 5 across all four canonical H-classes (Heisenberg, Ising, XY, soft XY+YX). The N-scan with γ_l = 0.05·(l+1) and the cross-class sanity at N = 4 with γ = [0.05, 0.10, 0.15, 0.20] both produced gap = candidate − truth = 0 to machine precision (~10⁻¹⁴). The verification ruled out the three failure hypotheses surfaced in the Phase 1 planning document:
- **(i) N=3-only artifact.** Refuted: formula is bit-exact at N = 3, 4, 5; no overlap-bond defect.
- **(ii) general-formula gap at all N.** Refuted: gap is 0 to machine precision.
- **(iii) centering mismatch.** Refuted: F1-centered L_Dc = L_D + σ · I matches the architect's assumed centering, and is the centering that gives the closed form below.

Phase 2 promotes the formula to a typed Tier-1-derived claim [`F49NonUniformCrossTermClaim`](../../compute/RCPsiSquared.Core/F1/F49NonUniformCrossTermClaim.cs) on the `F1KnowledgeBase` and closes the "Open follow-ups" item in [PROOF_F1_NONUNIFORM_GAMMA](PROOF_F1_NONUNIFORM_GAMMA.md).

## Proof

### Step 1: Anti-commutator decomposition via diagonal L_Dc

L_Dc is diagonal in the Pauli-string basis: `L_Dc · σ_α = d_α^c · σ_α` with `d_α^c = Σ_l γ_l · ε_l(α)` (Conventions). Hence the anti-commutator with the Hamiltonian superoperator inherits a pointwise structure:

    {L_H, L_Dc}_{αβ}  =  (L_H · L_Dc + L_Dc · L_H)_{αβ}
                       =  (L_H)_{αβ} · d_β^c  +  d_α^c · (L_H)_{αβ}
                       =  (L_H)_{αβ} · (d_α^c + d_β^c).

(Same as the F49 uniform proof, [Step 5 of PROOF_CROSS_TERM_FORMULA](PROOF_CROSS_TERM_FORMULA.md), with `d_α^c + d_β^c` replacing `2γ · (N − w_a − w_b)`.) Therefore

    ‖{L_H, L_Dc}‖²_F  =  Σ_{αβ} |(L_H)_{αβ}|² · (d_α^c + d_β^c)².

### Step 2: Bond-site separation under a single bond's transition

Restrict to a single bond term `H_b` acting on bond sites (i, j) := (i_b, j_b). By [PROOF_CROSS_TERM_FORMULA Lemma 3](PROOF_CROSS_TERM_FORMULA.md), every nonzero Pauli-basis transition of `L_H^bond_b` changes both bond sites and leaves all spectator sites unchanged: α_m = β_m for every m ∉ {i, j}. The per-site signs ε_l decompose accordingly:

    d_α^c + d_β^c
       =  Σ_l γ_l · ε_l(α)  +  Σ_l γ_l · ε_l(β)
       =  Σ_{m ∉ {i, j}} γ_m · (ε_m(α) + ε_m(β))            (spectators agree, ε_m(α) = ε_m(β))
        + γ_i · (ε_i(α) + ε_i(β))                            (bond site i)
        + γ_j · (ε_j(α) + ε_j(β))                            (bond site j)
       =  2 · Σ_{m ∉ {i, j}} γ_m · ε_m(α)
        + γ_i · A  +  γ_j · B,

where A := ε_i(α) + ε_i(β) and B := ε_j(α) + ε_j(β). Each of A, B ∈ {−2, 0, +2}.

**Bond-sum rule (Π²-class refinement of F49 Lemma 2):** For every nonzero (L_H^bond_b)_{αβ}, A + B is fixed by the bond Pauli class:
- ZZ-class bonds (both bond Paulis in {X, Y} or both in {I, Z}, the F49-shadow-balanced class): A + B = 0. This is the bond-sum rule of [PROOF_CROSS_TERM_FORMULA Lemma 2](PROOF_CROSS_TERM_FORMULA.md), now read at the ε-signed level instead of the w_XY-level: w_XY(α) + w_XY(β) = 2 on the bond corresponds to ε_i(α) + ε_i(β) + ε_j(α) + ε_j(β) = (1 − 2w_XY^{ij}(α)) + (1 − 2w_XY^{ij}(β)) = 2 − 2·2 = −2 + 2 = 0. (Equivalent statements.)
- Shadow-crossing bonds (e.g., XZ, YZ) violate the bond-sum rule; this is the F49c case and is outside the scope of the present proof.

For shadow-balanced bonds (the only case we consider here), B = −A. Substituting:

    d_α^c + d_β^c  =  2 · Σ_{m ∉ {i, j}} γ_m · ε_m(α)  +  A · (γ_i − γ_j).

### Step 3: Squaring and spectator-bond independence

Squaring:

    (d_α^c + d_β^c)²
       =  4 · (Σ_{m ∉ {i, j}} γ_m · ε_m(α))²
        + A² · (γ_i − γ_j)²
        + 4 · A · (γ_i − γ_j) · Σ_{m ∉ {i, j}} γ_m · ε_m(α).

The cross-term vanishes after summation over α (fixing the bond letters). The spectator-site ε_m(α) values are independent of the bond-site letters (Lemma 3: bond transitions leave spectators unchanged), so for each fixed bond-Pauli pair the sum Σ_α (with fixed bond letters, varying spectator letters) factorises:

    Σ_α |(L_H^bond_b)_{αβ}|² · 4 · A · (γ_i − γ_j) · Σ_{m ∉ {i, j}} γ_m · ε_m(α)
       =  Σ_{bond letters} |...|²_{bond} · 4A(γ_i−γ_j) · Σ_{spectator letters} Σ_m γ_m · ε_m(α)
       =  0

because Σ_{spectator letters at site m} ε_m(α) = (+1) + (−1) + (−1) + (+1) = 0 (the four Pauli letters I, Z, X, Y at site m give signs +1, +1, −1, −1, summing to zero; same vanishing as the cross-term of [PROOF_CROSS_TERM_FORMULA Lemma 1](PROOF_CROSS_TERM_FORMULA.md)). Hence

    Σ_{αβ} |(L_H^bond_b)_{αβ}|² · (d_α^c + d_β^c)²
       =  4 · Σ_{αβ} |(L_H^bond_b)_{αβ}|² · (Σ_{m ∉ {i, j}} γ_m · ε_m(α))²
        + (γ_i − γ_j)² · Σ_{αβ} |(L_H^bond_b)_{αβ}|² · A².

### Step 4: Spectator-sum reduces to (N−2)-version of F49 with weighting Σγ_m²

Expand the spectator square:

    (Σ_{m ∉ {i, j}} γ_m · ε_m(α))²
       =  Σ_{m ∉ {i, j}} γ_m² · ε_m(α)²  +  Σ_{m ≠ m', both ∉ {i, j}} γ_m γ_{m'} · ε_m(α) ε_{m'}(α).

Average over spectator letters (which run independently over the 4 Pauli choices, by Lemma 3):
- ε_m(α)² = 1 always (diagonal sum) ⟹ Σ_α ε_m(α)² = 4^N over all strings (the 4 spectator letters contribute uniformly).
- For m ≠ m', the spectator letters at sites m, m' are independent: Σ_{α_m} ε_m(α) · Σ_{α_{m'}} ε_{m'}(α) = 0 · 0 = 0 (each per-site sum vanishes; same calculation as F49 Lemma 1 / PROOF_CROSS_TERM_FORMULA Step 1).

The cross-spectator term vanishes by independence. The diagonal-spectator term, weighted by |(L_H^bond_b)_{αβ}|² (which depends only on the bond Paulis and the bond-target letters, fixed under spectator variation), reduces to

    Σ_{αβ} |(L_H^bond_b)_{αβ}|² · (Σ_{m ∉ {i, j}} γ_m · ε_m(α))²
       =  (Σ_{m ∉ {i, j}} γ_m²) · Σ_{αβ} |(L_H^bond_b)_{αβ}|²
       =  ‖L_H^bond_b‖²_F · Σ_{m ∉ {i, j}} γ_m².

(This is the per-site weighted analogue of [PROOF_CROSS_TERM_FORMULA Step 3 "spectator variance = N−2"](PROOF_CROSS_TERM_FORMULA.md). The uniform-γ identity `Σ_m γ² = (N−2)γ²` over spectators is replaced by `Σ_{m ∉ {i, j}} γ_m²`.)

Therefore the spectator part of the per-bond cross-term contribution is

    4 · ‖L_H^bond_b‖²_F · Σ_{m ∉ {i, j}} γ_m².

### Step 5: A-classification of bond Pauli pairs gives G(bond_b, H)

Define G(bond_b, H) := Σ_{αβ} |(L_H^bond_b)_{αβ}|² · A². The asymmetry contribution at this bond is

    G(bond_b, H) · (γ_i − γ_j)².

To compute G, decompose the bond Pauli pairs by the value of A ∈ {−2, 0, +2}. For a single bond term `α_i ⊗ β_j` (single Pauli pair, no sum), the commutator action `[α_i ⊗ β_j, P_i ⊗ Q_j]` on bond Paulis P_i, Q_j produces a transition with definite ε-shift A at site i: `A = ε_i(target) + ε_i(source) = ε_i(α · P) + ε_i(P) = ε_i(P) · (ε_i(α · P) / ε_i(P) + 1)`. Tabulating: the bond commutator changes the bond letters according to the Pauli structure constants, and one can enumerate the 16 (P_i, Q_j) configurations directly.

**Key A-classification (proven by direct enumeration over the 16 Pauli-pair transitions at a single bond):** for a bond term `α_i β_j` with α, β ∈ {X, Y, Z}, every nonzero Pauli-basis transition has

    A  =  ±2   if  (α, β) ∈ ZZ-class (both in {I, Z}, i.e., ZZ for genuine bond terms)
    A  =  0    if  (α, β) ∈ XY-class (XX, YY, XY, YX terms; all four bond-Paulis-in-{X, Y} configurations)

The structural reason: a ZZ bond term is the commutator generator that flips bit_a on both bond sites simultaneously (since Z commutes with I, Z and anti-commutes with X, Y, so its commutator action takes Q_j ∈ {X, Y} ↔ Q_j' ∈ {Y, X} respectively, both ⊂ {X, Y}, keeping bond-site bit_a class). The non-trivial transition `(P_i Q_j) → (α P_i, β Q_j)` for ZZ has α_i acting as Z_i (flipping the bond bit_a of P_i only inside {X, Y} ↔ {X, Y} or {I, Z} ↔ {I, Z}, leaving ε_i invariant on P_i = I or Z, flipping it on P_i = X or Y), but with both ε-shifts at sites i, j carrying the SAME sign (`ε_i(α P_i) − ε_i(P_i)` and `ε_j(β Q_j) − ε_j(Q_j)` have the same sign because Z preserves bit_a-class of P_i, Q_j ∈ {X, Y} together or {I, Z} together via the bond-sum rule). Hence A = ε_i(source) + ε_i(target) attains ±2 (both signs at the same letter).

For an XY-class bond term (XX, YY, XY, YX, i.e., both α, β ∈ {X, Y}), the commutator action flips bit_a on the affected bond site: XX takes (I_i, I_j) → (X_i, X_j) (no, XX in commutator: `[X_i X_j, P_i Q_j] = X_i [X_j, Q_j] · P_i + [X_i, P_i] · X_j Q_j`). On bond Paulis (P_i, Q_j) ∈ {I, X, Y, Z}², `[X_i, P_i]` is nonzero iff P_i ∈ {Y, Z} (anti-commuting), changing P_i ∈ {Y, Z} ↔ {Z, Y}. The output P_i' lies in the OPPOSITE bit_a-class from the input: X·Y = iZ has bit_a = 0 (Z), X·Z = −iY has bit_a = 1 (Y). So XX flips bit_a on both sites, giving A = ε_i(source) + ε_i(target) = 0 (since ε flips, one is +1 and the other is −1). The same holds for YY, XY, YX. Hence A = 0 for all XY-class bond terms.

Combining: only the ZZ-fraction of `H_b` contributes to G(bond_b, H). For a bond Hamiltonian `H_b = Σ_t c_t · α_t^{(i)} β_t^{(j)}` decomposed into ZZ-class and XY-class parts:

    G(bond_b, H)  =  4 · ‖L_{(ZZ-class part of H_b)}^bond_b‖²_F.

(The factor 4 comes from A² = 4 on every ZZ-class transition; the XY-class transitions contribute A² = 0.)

For canonical H-classes:
- **Heisenberg J·(XX + YY + ZZ).** The bond decomposes as XX + YY (XY-class, 2/3 of the bond norm) + ZZ (1/3 of the bond norm). Hence `‖L_{ZZ-part}^bond‖² = (1/3) · ‖L_H^bond‖²` and `G(bond, H) = 4·(1/3)·‖L_H^bond‖² = (4/3) · ‖L_H^bond‖²`.
- **Ising J·ZZ.** The bond is entirely ZZ-class, so `‖L_{ZZ-part}^bond‖² = ‖L_H^bond‖²` and `G(bond, H) = 4 · ‖L_H^bond‖²`.
- **XY J·(XX + YY).** No ZZ component, so `G(bond, H) = 0`.
- **Soft Π²-odd J·(XY + YX).** Both terms are XY-class (one bond letter in {X, Y}, other in {X, Y}), so `G(bond, H) = 0`.

### Step 6: Multi-bond sum via disjoint-bond-supports lemma

By [PROOF_CROSS_TERM_FORMULA Lemma 3 Corollary](PROOF_CROSS_TERM_FORMULA.md) ("disjoint bond supports"), for any two bonds e, e' (including overlapping bonds that share a site), the Pauli-basis transition supports of L_H^e and L_H^{e'} are disjoint: no (α, β) pair receives nonzero contributions from both bonds. Hence

    ‖L_H‖²_F   =  Σ_b ‖L_H^bond_b‖²_F
    Σ_{αβ} |(L_H)_{αβ}|² · (d_α^c + d_β^c)²  =  Σ_b Σ_{αβ} |(L_H^bond_b)_{αβ}|² · (d_α^c + d_β^c)².

Combining Steps 3, 4, 5 inside each per-bond sum, then summing over bonds:

    ‖{L_H, L_Dc}‖²_F
       =  4 · Σ_b ‖L_H^bond_b‖²_F · Σ_{m ∉ bond_b} γ_m²
        + Σ_b G(bond_b, H) · (γ_{i_b} − γ_{j_b})².

This is the theorem.    ∎

**Recovery of uniform-γ formula.** Set γ_l ≡ γ for all l. The bond-asymmetry part vanishes (γ_i − γ_j = 0 on every bond). The spectator part collapses to

    Σ_b ‖L_H^bond_b‖²_F · Σ_{m ∉ bond_b} γ_m²
       =  γ² · Σ_b ‖L_H^bond_b‖²_F · (N − 2)
       =  γ² · (N − 2) · ‖L_H‖²_F

(using |spectators per bond| = N − 2 and disjoint bond supports). Multiplied by the prefactor 4 in the theorem, this gives the F49 uniform identity `4γ²·(N−2)·‖L_H‖²_F`. ✓

## Verification

[`simulations/f49_nonuniform_gamma_crossterm_verify.py`](../../simulations/f49_nonuniform_gamma_crossterm_verify.py) verifies the closed form in four sections plus a final assertion block:

1. **Multi-centering at N = 3 Heisenberg γ = [0.1, 0.2, 0.3].** Four candidate L_Dc centerings; identifies F1-centered `L_D + σ·I` as the centering that gives the theorem's value 163.84.
2. **Architect-candidate formula at the F1 centering.** Reports per-bond spectator + asymmetry breakdown; the closed form predicts 163.84 (bit-exact).
3. **N-scan, Heisenberg chain, γ_l = 0.05·(l+1) at N = 3, 4, 5.** All gaps `candidate − truth = 0` to machine precision: 40.96, 737.28, 8437.76 in turn.
4. **Cross-H-class sanity at N = 4 γ = [0.05, 0.10, 0.15, 0.20].** All four classes match bit-exact:
   - Heisenberg: 737.28 (G factor 4/3 on bond norm 1536)
   - Ising ZZ: 256.00 (G factor 4 on bond norm 512)
   - XY (XX+YY): 481.28 (G factor 0, spectator-only)
   - Soft XY+YX: 481.28 (G factor 0, spectator-only)
5. **Final assertion block.** All N ∈ {3, 4, 5} × {Heisenberg, Ising, XY, soft XY+YX} candidate-vs-truth pairs assert within absolute tolerance 1e-10. Script exits 0 on success; AssertionError surfaces any future regression.

All verifications pass at machine precision. Summary of per-bond Frobenius norms (J = 1) used by the closed form:

| N | Heisenberg | Ising ZZ | XY (XX+YY) | XY+YX soft |
|---|-----------:|---------:|-----------:|-----------:|
| 3 |        384 |      128 |        256 |        256 |
| 4 |       1536 |      512 |       1024 |       1024 |
| 5 |       6144 |     2048 |       4096 |       4096 |

(Each step multiplies by 4 because one additional spectator site contributes `tr(I_4) = 4` to the Frobenius-norm tensor calculation; intrinsic local 2-qubit bond superoperator norms are 96 / 32 / 64 / 64 respectively.)

## Diagnostic interpretation

The closed form makes `‖{L_H, L_Dc}‖²_F` a quantitative two-piece diagnostic under non-uniform γ:

- **Spectator part dominates at large N.** For uniform-γ scaling at fixed H and γ̄, the spectator part scales as `4γ̄² · (N−2) · ‖L_H‖²` (recovering the F49 uniform identity); the asymmetry part is `Σ_b G(bond_b, H) · (γ_i − γ_j)²` and is independent of N for fixed bond structure. The ratio spectator/asymmetry grows linearly in N.
- **Bond-asymmetry is the H-class fingerprint.** `G(bond_b, H) / ‖L_H^bond_b‖²` reads out the ZZ-fraction of the bond Hamiltonian: 4/3 for Heisenberg (ZZ is 1/3), 4 for Ising (ZZ is 100%), 0 for XY-class (no ZZ). Non-uniform γ extraction at fixed H can be inverted via the asymmetry signature: measure `‖{L_H, L_Dc}‖²` at multiple γ patterns, fit the (γ_i − γ_j)² coefficients per bond, recover the per-bond ZZ-content.
- **Bit_a-class preservation is the structural source.** The bond-asymmetry coefficient G(bond_b, H) = 4·‖L_{ZZ-part}^bond_b‖² traces to a simple structural fact: ZZ bond terms preserve bit_a-class at each site (Z keeps {I, Z} ⊂ {I, Z} and {X, Y} ⊂ {X, Y}), so the per-site ε-signs add coherently (A = ±2). XY-class terms flip bit_a at both sites simultaneously, so the per-site ε-signs cancel (A = 0). The bond-asymmetry sensitivity to (γ_i − γ_j) is therefore a direct readout of how much of `H_b` is bit_a-class-preserving.
- **Uniform-γ recovers F49.** At γ_l ≡ γ, the formula collapses to `4γ² · (N−2) · ‖L_H‖²_F`. The asymmetry vanishes, and the spectator part absorbs the entire F49 numerator. Conversely, the F49 cross-term ratio `R(N)² = (N−2)/(N·4^(N−1))` is preserved at uniform γ; for non-uniform γ, R depends on the specific γ pattern through both `Σ_b ‖L_H^bond_b‖² · Σ_{m ∉ bond_b} γ_m²` and `Σ_b G(bond_b, H) · (γ_i − γ_j)²`, so a single dimensionless R no longer exists.

## Cross-references

### Repository entries

- **F1 palindrome identity** ([`docs/ANALYTICAL_FORMULAS.md` F1](../ANALYTICAL_FORMULAS.md#f1-palindrome-equation-tier-1-proven), [MIRROR_SYMMETRY_PROOF.md](MIRROR_SYMMETRY_PROOF.md)): the parent identity Π·L·Π⁻¹ + L + 2σ·I = 0 whose σ-shift gives L_Dc.
- **F49 uniform cross-term formula** ([`docs/ANALYTICAL_FORMULAS.md` F49](../ANALYTICAL_FORMULAS.md#f49-cross-term-formula-tier-1-proven), [PROOF_CROSS_TERM_FORMULA.md](PROOF_CROSS_TERM_FORMULA.md)): the uniform-γ parent identity. The present proof is its non-uniform γ extension; Steps 1, 2, 3, 4, 6 use the same bond-sum rule (Lemma 2), spectator-variance reduction (Step 3), and disjoint-bond-supports lemma (Lemma 3 + Corollary).
- **F1 non-uniform γ H-block closure** ([PROOF_F1_NONUNIFORM_GAMMA.md](PROOF_F1_NONUNIFORM_GAMMA.md)): sibling closure showing the H-block residual scaling factor F(N, G) is γ-independent. The "Open follow-ups" of that proof noted the F49 cross-term gap at non-uniform γ; the present proof closes that follow-up.
- **F1 T1 block closed form** ([PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md](PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md), [F1T1ResidualClosedForm](../../compute/RCPsiSquared.Core/F1/F1T1ResidualClosedForm.cs)): per-Pauli-class dissipator-block closed form; H-block-independent. Parallel to the present proof in that both extract a per-bond / per-site analytic structure from a Pauli-basis tensor calculation.
- **F1 depol block closed form** ([PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md](PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md), [F1DepolResidualClosedForm](../../compute/RCPsiSquared.Core/F1/F1DepolResidualClosedForm.cs)): another sibling dissipator-block closed form.
- **F49c shadow-crossing companion** ([`docs/ANALYTICAL_FORMULAS.md` F49c](../ANALYTICAL_FORMULAS.md), [PROOF_CROSS_TERM_CROSSING.md](PROOF_CROSS_TERM_CROSSING.md)): companion formula for the case where bond Pauli pairs violate the bond-sum rule. Non-uniform extension to shadow-crossing couplings is OUT OF SCOPE for the present proof; the theorem above applies only to shadow-balanced couplings (Heisenberg, Ising, XY family).

### Typed claims

- [`compute/RCPsiSquared.Core/F1/F49NonUniformCrossTermClaim.cs`](../../compute/RCPsiSquared.Core/F1/F49NonUniformCrossTermClaim.cs): Tier-1-derived typed claim for this closed form, registered on `F1KnowledgeBase`. Closes the "Open follow-ups" item of [PROOF_F1_NONUNIFORM_GAMMA](PROOF_F1_NONUNIFORM_GAMMA.md).
- [`compute/RCPsiSquared.Core/F1/PalindromeResidualScalingClaim.cs`](../../compute/RCPsiSquared.Core/F1/PalindromeResidualScalingClaim.cs): companion H-block scaling claim, γ-independent (proved by [PROOF_F1_NONUNIFORM_GAMMA](PROOF_F1_NONUNIFORM_GAMMA.md)).

### Scripts

- [`simulations/f49_nonuniform_gamma_crossterm_verify.py`](../../simulations/f49_nonuniform_gamma_crossterm_verify.py): verification script; Phase 1 exploratory + Phase 2 assertion block.

### Memory

- `project_palindrome_frobenius_scaling`: the F1 H-block scaling memory; non-uniform γ extension for the cross-term closed here.
