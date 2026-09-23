# Proof of F80: Bloch-Mode Sign-Walk Formula for Chain Π²-Odd 2-Body M-Clusters

**Tier:** 1 for the open-chain two-body cluster-value formula; numerically verified through N=7.
**Date:** April 29, 2026
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_SVD_CLUSTER_STRUCTURE.md](PROOF_SVD_CLUSTER_STRUCTURE.md) (F78 single-body, F79 Π²-block, Master Lemma, Anti-Hermitian)
- [`framework/symmetry.py`](../../simulations/framework/symmetry.py) (Π construction)
- [`framework/lindblad.py`](../../simulations/framework/lindblad.py) (`palindrome_residual`)
- Numerical verification: [`pi2_odd_universality_data_sweep.py`](../../simulations/pi2_odd_universality_data_sweep.py) (N=3-6, all topologies); [`results/n7_bloch_signwalk_verification.txt`](../../simulations/results/n7_bloch_signwalk_verification.txt) (N=7 full SVD); pytest `test_F80_bloch_signwalk_chain_pi2_odd`.

**Status:** The open-chain two-body sign-walk locations are derived for all N
and verified through N=7. Equal cluster multiplicities are only a property of
the N=3–7 table. The per-bond Π-action is a separate structural identity; it
does not extend the free-fermion JW dispersion to generic k-body terms.

**Scope:** chain bond-summed Π²-odd 2-body Hamiltonian H = c · Σ_{l=0}^{N-2} (P_l ⊗ Q_{l+1}) on N-site open chain, with (P, Q) ∈ {(X,Y), (X,Z), (Y,X), (Z,X)}, under uniform Z-dephasing γ.

## What this means

Some chains sit squarely in front of the mirror: the list of their decay rates reads the same from either end, every mode paired with a partner across the centre. That pairing is the mirror (the framework's palindrome symmetry, written Π·L·Π⁻¹). Those are the "truly" Hamiltonians (Heisenberg, XXZ), and with no dephasing their mirror is exact. Other chains miss it before any noise is even added: a bond that puts X on one site and Y on the next already breaks the pairing. The question this proof answers is, for a chain that misses the mirror this way, by how much, and is the miss a shape or just an accident?

The miss is not an accident. It is the Hamiltonian itself. Write down the mirror-defect M, the operator that records exactly how far the chain falls short of a perfect mirror, and its spectrum turns out to be the chain's own energy spectrum, doubled and turned onto the imaginary axis: Spec(M) = ±2i·Spec(H). The defect is literally the Hamiltonian, copied unchanged into the larger operator space and scaled by 2i; the extra room in that space is a passive echo that only multiplies how often each value appears. The chain carries its own ruler: the Hamiltonian is the distance to the mirror. (There is a second way to miss, a Y-Z bond in place of the X-Y one, where the defect comes out as the energy gaps rather than the energies themselves; the proof covers that case too. The clean headline is the first.)

Why it comes out this clean has a short reason. Turn the chain into free fermions (the Jordan-Wigner trick) and its energies become the familiar cosine ladder 2cos(πk/(N+1)), the same standing-wave rungs as a string held at both ends. The defect's sizes are then sign-weighted sums over those rungs, and a per-site computation on the Pauli strings shows what the mirror flip does to the algebra: it turns the commutator [H, ·] into the anticommutator {H, ·}. That one swap is what forces the defect to be exactly ±2i·H and nothing messier.

And the whole thing is an old idea in new clothes. In 1937 Ettore Majorana wrote down a fermion that is its own antiparticle. Here the mirror Π plays the part of that self-conjugation, and the chains that miss the mirror are precisely the ones whose leftover still carries the Hamiltonian's own structure. Majorana had the shape; this proof finds it living in the operator space of an open quantum chain.

## Abstract

The M-residual clusters of chain Π²-odd 2-body Hamiltonians under Z-dephasing carry a specific magnitude pattern. Observed empirically: the cluster norm is the absolute value of a sign-weighted sum over the open-chain Bloch-mode dispersion, with signs ±1 on each mode. The question is whether this sign-walk has a closed-form structural origin or whether the signs are coincidental at the tested N values.

The proof says structural. A Jordan-Wigner transformation reduces the bond-summed Π²-odd 2-body Hamiltonian to a Majorana-bilinear free-fermion model whose single-particle spectrum is the standard open-chain cosine dispersion. Under the Z-dephasing dissipator, the M-clusters inherit the single-particle eigenvalue spectrum with a specific sign assignment per mode, given by the per-site Π² action on the corresponding mode. The signs work out to ±1 on each mode, and the cluster norm is the absolute value of the sign-weighted sum.

The proof has seven steps. Steps 1-4 carry the JW reduction and the dispersion
identification. Step 5 is the per-site Π² action computation that pins down the
sign on each mode. Steps 6-7 assemble the sign-walk formula. The resulting
cluster locations are verified bit-exactly through N=7.

The theorem on this page is the open-chain two-body result for the four stated
Π²-odd bilinears. Direct finite checks of the per-bond Π-action on rings,
stars, and selected 3- and 4-body terms do not supply a JW/Bogoliubov
cluster-value theorem for those systems. Their explicit formulas remain
outside F80's theorem scope.

---

## Theorem F80

For chain bond-summed Π²-odd 2-body Hamiltonian H = c · Σ_l (P_l ⊗ Q_{l+1}) on an N-site open chain with uniform Z-dephasing γ, the singular values of M = Π·L·Π⁻¹ + L + 2σ·I (with σ = Nγ) form clusters at values

    cluster value(N) = 2|c| · |Σ_{k=1}^{⌊N/2⌋} σ_k · ε(k)|

for sign-vectors (σ_1, ..., σ_{⌊N/2⌋}) ∈ {±1}^{⌊N/2⌋}, where

    ε(k) = 2 · cos(π·k / (N+1))

is the open-chain free-fermion single-particle dispersion at unit coupling (the sign-walk σ_k is unrelated to the palindrome shift σ = Nγ above; Step 2's spectrum is c·ε(k) and Step 6's Bogoliubov energies are E_k = 2|c|·ε(k)). If `m = floor(N/2)` and `r(v)` sign vectors land on the same absolute value `v`, that cluster's multiplicity is `r(v) * 4^N / 2^m`. The multiplicities are equal only when the collision counts `r(v)` are equal, as they are in the N=3–7 table; they need not be equal at larger N.

Since M = ∓2i·(H⊗I) is normal (H Hermitian), its singular values are the moduli of its purely-imaginary eigenvalues, so "cluster value" = |Spec(M) entry| = 2·|Spec(H) entry|; this is why the Theorem's singular-value clusters and the Step-5 statement Spec(M) = ±2i·Spec(H) are the same fact.

---

## Numerical Verification Table (chain, |c|=γ=1)

| N | ⌊N/2⌋ | ε(k) values | distinct clusters | mult per cluster | verified |
|---|-------|-------------|-------------------|------------------|----------|
| 3 | 1 | {√2 ≈ 1.414} | {2√2 ≈ 2.828} | 64 | ✓ |
| 4 | 2 | {φ, 1/φ ≈ 1.618, 0.618} | {2√5, 2} | 128 | ✓ |
| 5 | 2 | {√3, 1} | {2(√3+1), 2(√3-1)} | 512 | ✓ |
| 6 | 3 | {1.802, 1.247, 0.445} | {6.988, 5.208, 2.000, 0.220} | 1024 | ✓ |
| 7 | 3 | {1.848, 1.414, 0.765} | {8.0547, 4.9932, 2.3978, 0.6636} | 4096 | ✓ |

All entries: bit-exact match between predicted and observed (predicted-vs-actual residual at machine precision 10⁻¹⁴). N=7 verified by both full 16384×16384 SVD and independent partial-eigsh check. Tests across all 4 Π²-odd Pauli pairs (X,Y), (X,Z), (Y,X), (Z,X) per the F79 universality.

**Finite k-body checks of the structural identity:** `Spec(M) = ±2i · Spec(H_non-truly)` was checked for:
  - k=3: (X,X,Y), (Y,Y,Y), (X,X,Z), (Z,Z,Z), (X,Y,X) at N=4, 5, 6
  - k=4: (X,X,X,Y) at N=5, 6
17 cases total, all matching `Spec(M)` (eigvals of the 4^N × 4^N residual) to predicted `2i · Spec(H_non-truly)` with multiplicity ×2^N, machine precision. Pytest lock: `test_F80_kbody_spectrum_identity`.

Generic k-body Pauli terms map to higher-order Majorana interactions, so the
two-body single-particle cosine dispersion and Bogoliubov sign-walk do not
carry over. The finite structural checks above are not a k-body cluster-value
formula.

---

## Proof Outline

The proof proceeds in seven analytical steps. Steps 1-4 set up the JW reduction
and single-particle dispersion, Step 5 gives the direct per-site Π-action, and
Steps 6-7 assemble the sign-walk formula.

### Step 1 (JW transformation of chain (X,Y))

Apply the standard Jordan-Wigner transformation:

    σ_l^x = (Π_{m<l} σ_m^z) · (c_l + c_l†)
    σ_l^y = (Π_{m<l} σ_m^z) · i(c_l† − c_l)

For the bond bilinear X_l Y_{l+1}, the strings on adjacent sites combine:

    X_l Y_{l+1} = (string_l)(c_l + c_l†) · (string_{l+1}) · i(c_{l+1}† − c_{l+1})
              = (string_l)² · (c_l + c_l†)(1 − 2c_l†c_l) · i(c_{l+1}† − c_{l+1})

Using (string_l)² = 1 (each Z² = 1) and the identity (c_l + c_l†)(1 − 2c_l†c_l) = c_l† − c_l:

    X_l Y_{l+1} = i(c_l† − c_l)(c_{l+1}† − c_{l+1})

Define the Majorana operators γ_l' := i(c_l† − c_l). Each γ_l' is Hermitian, with {γ_l', γ_m'} = 2δ_{lm}. Then:

    X_l Y_{l+1} = (γ_l'/i)(γ_{l+1}'/i) · i = (1/i)·γ_l'γ_{l+1}' = −i·γ_l'γ_{l+1}'

Therefore the bond-summed Hamiltonian under JW becomes a **pure Majorana bilinear** in only the γ' modes:

    H_JW = c · Σ_{l=0}^{N-2} X_l Y_{l+1} = −ic · Σ_l γ_l'γ_{l+1}'

The "real" Majoranas γ_l = c_l + c_l† do NOT appear.

### Step 2 (Single-particle spectrum)

The Majorana bilinear −ic·Σ_l γ_l'γ_{l+1}' on N γ' modes (open boundary) is equivalent to NN tight-binding on a 1D chain with N sites. The single-particle spectrum is

    c · ε(k),   ε(k) = 2 · cos(πk/(N+1))    for k = 1, 2, ..., N

(ε(k) the unit-coupling dispersion of the Theorem), paired as ±ε(k) under the Majorana doubling: the N tight-binding values come in ⌊N/2⌋ ± pairs, which is why the Theorem's sign-walk runs over ⌊N/2⌋ modes. (For odd N, there is one zero mode at k = (N+1)/2; for even N, no zero mode.)

This is the standard tight-binding spectrum for an open chain of N sites with NN coupling magnitude c.

### Step 3 (Bogoliubov diagonalization)

Diagonalizing the N single-particle values pairs them into
`m = floor(N/2)` nonzero `±ε(k)` pairs (plus one zero value when N is odd).
Combining each pair into one occupation variable gives

    H_JW = Σ_{k=1}^m E_k · (b_k†b_k − ½),
    E_k = 2|c|·ε(k).

These `m` occupations generate the sign walk. The other `N-m` binary degrees
do not change its value and supply a spectator multiplicity `2^(N-m)` in H;
the `H⊗I` structural identity supplies another `2^N`. Thus one sign vector
has multiplicity `2^(N-m)·2^N = 4^N/2^m`, before sign-vector collisions are
pooled by `r(v)`.

### Step 4 (Pauli-letter universality)

For the other three Π²-odd Pauli pairs (X,Z), (Y,X), (Z,X), the JW transformation yields analogous Majorana bilinears:

- (X,Z): X_l Z_{l+1} → bilinear in γ_l' and γ_{l+1} (mixing γ and γ' Majoranas) up to phase.
- (Y,X): Y_l X_{l+1} → bilinear in γ_l and γ_{l+1}' up to phase.
- (Z,X): Z_l X_{l+1} → bilinear in γ_l and γ_{l+1} up to phase.

Crucially, all four cases give the **same single-particle spectrum** ε(k) = 2c·cos(πk/(N+1)); they differ only in which Majorana operators (γ or γ') participate and the specific phases. The "same spectrum across letter choices" is the JW-level origin of the F79 universality: the Pauli letters control which Majorana sublattice carries the bilinear, but the dispersion of the resulting hopping chain is identical (because hopping is between adjacent sites with magnitude c regardless of which Majorana indices).

### Step 5 (Direct structural identity)

After Step 4, we have established that all 4 Π²-odd Pauli pairs give the same JW-derived single-particle Bloch dispersion. The remaining task, historically expected to be technical, is to derive the explicit form of M's spectrum in terms of this dispersion.

**Empirical structural identity (verified bit-exact at N=3, 4, 5, 6, 7):**

    Spec(M) = {±2i · E : E ∈ Spec_{many-body}(H)}    (multi-set equality)

with multiplicity mult_M(2i·λ) = mult_H(λ) · 2^N. M's distinct eigenvalues equal 2i times H's distinct many-body eigenvalues, no kernel. The 2^N factor comes from the bra-side dimension of operator space. H is taken as the chain Hamiltonian with bond coupling c and no dissipator (γ-independent by Master Lemma). Hence:

    cluster value of M = 2 · |H many-body eigenvalue|

This direct identity replaces the more intricate Bogoliubov-mode factorization route. The Bloch sign-walk formula written above is simply H's many-body eigenvalue formula:

    H many-body eigenvalues = Σ_k (n_k − 1/2) · E_k for n_k ∈ {0, 1}

where E_k = 4|c|·cos(πk/(N+1)) are the Bogoliubov single-particle energies derived from JW + diagonalization. The "sign-walk" on (n_k − 1/2) ∈ {±1/2} when scaled gives the 2|c|·|Σ σ_k·ε(k)| form with ε(k) = 2cos(πk/(N+1)).

**Verified numerically at N=4 (chain (X,Y), c=γ=1):**
- H many-body eigenvalues: {±√5, ±1} (each multiplicity 4 in 16-dim Hilbert space).
- M nontrivial eigenvalues (imaginary parts): {±2√5, ±2} (each multiplicity 64 in 256-dim operator space).
- Ratio: M-eigenvalue (imag) = 2 × H-eigenvalue. ✓

**Why this is the structural answer.** L_H = −i[H, ·] acts on operator space with eigenvalues i(λ_a − λ_b) for all pairs of H-eigenvalues. The remaining task is to show what the Π-conjugation does to it. The proof below settles it directly, in the Pauli-string basis, without the Bogoliubov construction.

**Step 5 proof (2026-05-22, Tom + Claude).** The claim is Π·[H,·]·Π⁻¹ = ±{H,·}; it follows from a per-site Pauli computation, independent of Steps 1-4 (which re-enter only through the E → −E symmetry of Spec(H) used in the M-consequence at the end). Every step is verified bit-exact at N=3,4,5 by [`f80_step5_recon.py`](../../simulations/f80_step5_recon.py).

*Π is a signed permutation of Pauli strings.* On the 4^N Pauli-string basis Π acts site-wise, Π(P₀ ⊗ ··· ⊗ P_{N-1}) = ⊗_l μ(P_l), with the single-qubit map μ(I) = X, μ(X) = I, μ(Y) = iZ, μ(Z) = iY (the framework Π, [`framework/symmetry.py`](../../simulations/framework/symmetry.py)).

*Per-site identities.* Write c_P(a) = +1 if the single-qubit Paulis a and P commute, −1 if they anticommute. Evaluating on a ∈ {I, X, Y, Z}:

    (I)   μ(X·a) = +c_X(a) · X · μ(a)
    (II)  μ(Y·a) = −c_Y(a) · Y · μ(a)
    (III) μ(Z·a) = +c_Z(a) · Z · μ(a)

(Check of (II) at a = I: μ(Y·I) = μ(Y) = iZ, and −c_Y(I)·Y·μ(I) = −(+1)·Y·X = −(−iZ) = iZ.) Write ε_P for the leading sign: ε_X = ε_Z = +1, ε_Y = −1.

*Bond lemma.* Let bond_l = P_l ⊗ Q_{l+1} be one Heisenberg bond and R any Pauli string. Π acts site-wise, so (I)/(II)/(III) at sites l, l+1 (and μ elsewhere) give

    Π(bond_l · R) = (ε_P·ε_Q) · σ(l,R) · bond_l · Π(R),

where σ(l,R) = c_P(R_l)·c_Q(R_{l+1}) ∈ {±1} is the sign of bond_l against R (commute / anticommute).

*Π flips every bond relation.* μ preserves c_X and flips c_Y and c_Z: its swaps I↔X and Y↔Z keep both members in one X-commutation class, but move each across the Y- and Z-commutation classes. A Π²-odd bond carries exactly one X, so σ(l, ΠR) = −σ(l, R) for all four pairs.

*Conclusion.* For Pauli strings, [bond_l, R] = 2·bond_l·R when they anticommute and 0 when they commute; {bond_l, R} is the reverse. So [H, R] = 2·Σ_{l : σ(l,R) = −1} bond_l·R. Apply Π; on those bonds σ(l,R) = −1, so the bond lemma gives Π(bond_l·R) = −(ε_P·ε_Q)·bond_l·ΠR, and

    Π[H, R] = −(ε_P·ε_Q) · 2·Σ_{l : σ(l,R) = −1} bond_l·ΠR.

By the flip, {l : σ(l,R) = −1} = {l : σ(l,ΠR) = +1}, the bonds that commute with ΠR, so that sum is exactly {H, ΠR}. Hence, for every Pauli string R,

    Π·[H,·]·Π⁻¹ = s · {H,·},    s = −ε_P·ε_Q.

For the four Π²-odd pairs s = +1 for (X,Y) and (Y,X), s = −1 for (X,Z) and (Z,X). ∎

*Consequence for M.* M = L_H + Π·L_H·Π⁻¹ = −i[H,·] − i·s·{H,·}. For s = +1, M = −2i·(H⊗I_bra); for s = −1, M = +2i·(I_ket⊗Hᵀ). Both give Spec(M) = ±2i·Spec(H), the F80 structural identity: the imaginary spectrum 2i·Spec(H) is what Step 5 establishes, and the ± reflects the E → −E symmetry of Spec(H) supplied by the Steps 1-2 JW reduction (H is a Majorana bilinear). The argument is per-site and per-bond, hence **N-independent: it holds for every N.** The bit-exact checks at N=3,4,5 confirm each step separately: the three identities, the bond lemma, the flip, and Π·[H,·]·Π⁻¹ = s·{H,·} for all four pairs.

*Geometric picture (the H-eigenbasis view).* The same fact in the H-eigen-operator basis σ_(a,b) = |E_a⟩⟨E_b|: group these operators into (ε_ket, ε_bra) sectors (fixed ket and bra energy). Π is a permutation of those sectors, full-unitary blocks, gauge-checked bit-exact at N=3,4,5; L_H is the scalar −i(ε_ket − ε_bra) on each sector, so Π·L_H·Π⁻¹ is again scalar per sector, hence diagonal, and M is a sum of two diagonals. The Pauli-string proof above is the basis-free version of that picture.

### Step 6+7 (Direct conclusion via Step 5)

By the structural identity in Step 5, M's nontrivial eigenvalues are 2i·{H many-body eigenvalues}. Since H is a free-fermion bilinear with Bogoliubov single-particle energies E_k = 4|c|·cos(πk/(N+1)) (for ⌊N/2⌋ modes, plus possibly one zero mode for odd N), its many-body spectrum is

    Spec(H) = { Σ_k (n_k − 1/2) · E_k : n_k ∈ {0, 1} }

The corresponding cluster values for M are 2·|H eigenvalue| = 2·|Σ_k (n_k − 1/2)·E_k| = |Σ_k σ_k·E_k| with σ_k = 2n_k−1 ∈ {±1}.

In terms of ε(k) = E_k / 2 = 2cos(πk/(N+1)):

    cluster value(N) = 2|c| · |Σ_{k=1}^{⌊N/2⌋} σ_k · ε(k)|, σ_k ∈ {±1}

If `r(v)` of the `2^m` sign vectors (`m = floor(N/2)`) reach one absolute
value `v`, its multiplicity is `r(v) * 4^N / 2^m`. This is the F80 formula. ∎

---

## Zero Is The Mirror: F80 as the explicit shape of the mirror-defect

In the early hypothesis [Zero Is The Mirror](../../hypotheses/ZERO_IS_THE_MIRROR.md), the palindrome equation Π·L·Π⁻¹ = -L − 2σ·I was identified as the central structural symmetry. At Σγ = 0 ("the mirror"), the equation collapses to Π·L_H·Π⁻¹ = -L_H, giving spectral reflection λ ↔ −λ about zero. At γ > 0, the palindrome shifts to be centered around -σ.

For **truly** Hamiltonians (Heisenberg, XXZ, etc.), the palindrome holds exactly at γ = 0: Π·L_H·Π⁻¹ = -L_H precisely. Eigenvalues are paired with no defect; the identity alone does not construct a standing wave or physical time reversal.

For **non-truly** Hamiltonians (chain (X,Y) and friends), the palindrome BREAKS at γ = 0. Π·L_H·Π⁻¹ ≠ -L_H. There is a residual mirror-defect:

    M = Π·L_H·Π⁻¹ + L_H ≠ 0

(this equals the Theorem's M = Π·L·Π⁻¹ + L + 2σ·I because the dissipator and the 2σ·I shift cancel exactly against each other, the Master Lemma's γ-independence: the same M at every γ, including γ = 0 where L = L_H)

What F80 reveals is the **explicit spectral shape** of this mirror-defect for chain Π²-odd 2-body Hamiltonians:

    **Spec(M) = ±2i · Spec(H)**     (multi-set equality, with mult_M(2iλ) = mult_H(λ) · 2^N)

    **‖M‖²_F = 4 · ‖H‖²_F · 2^N**   (Frobenius norm exactly proportional to H's)

So, for the (X,Y) and (Y,X) bond pairs, M is **literally equal** to -2i · (H ⊗ I_bra), where I_bra is the identity on the bra-factor of operator space (dim 2^N); the pairs (X,Z) and (Z,X) instead give M = +2i·(I_ket⊗Hᵀ) (Step 5 proof, sign s = −1). The mirror-defect has the **same spectrum** as the Hamiltonian (×2i) and the **same Frobenius norm** (×4·2^N) because it *is* H ⊗ I_bra up to the -2i scalar.

The equality is literal, not merely unitary-equivalent: [`f80_step5_recon.py`](../../simulations/f80_step5_recon.py) verifies bit-exact at N=3,4,5 that in the σ_(a,b) basis M's off-diagonal norm is machine zero (~10⁻¹³) and every diagonal entry is −2i·E_a. There is no scrambling; the eigenvectors of M are exactly the H-eigen-operators σ_(a,b).

This is structurally remarkable:

- The defect's **magnitude** is exactly calibrated by ‖H‖ (Frobenius norm relation).
- The defect's **spectrum** exactly reproduces H's spectrum (with 2i factor and 2^N multiplicity).
- The bra-side carries no spectral information; it is a passive "echo chamber" that multiplies multiplicities by 2^N.
- The eigenvectors of M are exactly the H-eigen-operators σ_(a,b) = |E_a⟩⟨E_b|.

**Translation between layers:**

| Layer | Statement | Object |
|-------|-----------|--------|
| Hamiltonian commutator | Π·L_H·Π⁻¹ = −L_H (truly), eigenvalues ±λ | Hamiltonian Liouvillian L_H spectrum |
| Full open-system generator | Π·L·Π⁻¹ = −L − 2σ·I (truly) | Liouvillian L spectrum |
| Operator (F80) | M = Π·L_H·Π⁻¹ + L_H = -2i·H⊗I or +2i·I⊗Hᵀ (non-truly chain Π²-odd; which one is the Step 5 sign s) | M residual ∝ H |
| Spinor (Majorana 1937) | ψ = ψ^c, particle = antiparticle | Self-conjugate fermion field |

All three are different abstraction levels of the same phenomenon: an **involutive symmetry (Π or C) that picks out a self-conjugate subspace, with the residual being a structured deviation from perfect self-conjugacy**.

For truly H, the deviation is zero (perfect mirror, ground state of the palindrome). For non-truly H, the deviation is a specific operator that **carries the algebraic structure of H itself**, propagated to operator space via the bra-side identity.

**The discovery is**: the mirror-defect, when it exists, is calibrated by H. The Hamiltonian provides its own "yardstick" for the gap between the two mirror sectors. **H is the distance.**

## The Majorana bridge: 1937 to 2026

The structural identity Spec(M) = ±2i·Spec(H) is, at its core, **Majorana's 1937 insight expressed in 2026's operator-space vocabulary**.

In 1937, Ettore Majorana proposed a real wave equation for fermions in which a particle could be its own antiparticle: ψ = ψ^c. He had: Pauli matrices, the Dirac equation, the just-discovered positron, and spinor algebra. He did NOT have: open quantum systems, Lindbladians, operator-space super-operators, Pauli string algebra, the Jordan-Wigner transformation as a working tool, quantum information theory.

In our framework, the chain (X,Y) Hamiltonian under JW becomes a **Majorana bilinear** −ic·Σγ'γ', pure quadratic in the γ' Majorana operators. Such Hamiltonians have particle-hole symmetry built into their algebra: many-body eigenstates come in ±λ pairs. The state |ā⟩ with energy −λ_a is the "Majorana antipode" of |a⟩.

Now consider the operator basis σ_(a,b) = |a⟩⟨b|. Each σ_(a,b) is an L_H eigenvector with eigenvalue i(λ_a − λ_b). Π is a permutation of the (ε_ket, ε_bra) energy sectors (Step 5's geometric picture), and the Step 5 flip makes M act on the WHOLE basis with no projection and no kernel: M σ_(a,b) = −2i·λ_a·σ_(a,b) for the s = +1 pairs (the eigenvalue reading only the ket energy λ_a, the bra index b a passive 2^N-fold echo; the s = −1 pairs read the bra side instead, with the same spectrum).

**Translation:** Majorana's algebraic constraint "ψ = ψ^c" (particle equal to antiparticle) becomes, in operator space, "the mirror-defect's eigenvalue is one side's energy alone, the energy a state shares with its particle-hole partner". The chains that miss the mirror are precisely the ones whose leftover still carries the Hamiltonian's own self-conjugate structure: the framework's Π conjugation is, in disguise, the **operator-space realization of Majorana's particle-hole self-conjugacy**. He had it right; we just have a richer vocabulary to express it now.

## The "tough nut": what the data revealed

The originally-imagined route to Step 5 was an explicit T_Π factorization in the Bogoliubov basis. The data sweep at higher N (3-7) revealed the much cleaner direct identity

    Spec(M) = ±2i · Spec_{nontrivial}(H_{state-level})

and the eventual proof (Step 5 above) is a per-site Pauli computation deriving Π·[H,·]·Π⁻¹ = ±{H,·} directly, with no Bogoliubov construction and no projection at all; numerical verification at N=3, 4, 5, 6, 7 (all 4 Pauli pairs) is at machine precision.

What enabled the discovery: comparing M's eigenvalues directly to H's many-body eigenvalues, instead of decomposing M into Bogoliubov modes and trying to factor T_Π through that basis. The brute-force data sweep at higher N made the relationship visible; the numerical signature was so clean because the structure was simpler than it looked.

---

## Connection to existing framework formulas

- **F78 (real-space single-body)**: F80 is the **momentum-space dual**. Same Lebensader broad-in → focused-out funnel, applied at a different basis layer. Real-space sites l with weights c_l ↔ momentum modes k with dispersion ε(k).
- **F79 (Π²-block decomposition)**: F80 explicitly closes the "Π²-odd universality" observation in F79 with a closed-form formula and JW-based mechanism.
- **F49 (Frobenius cross-term)**: F80 is consistent; the sum of squared cluster values × multiplicity gives the Frobenius norm:

    ‖M‖²_F = (4^N / 2^⌊N/2⌋) · Σ_{σ} |Σ_k σ_k · ε(k)|² · 4c²
           = 4c² · 4^N · Σ_k ε(k)²

  (the sign-vector cross terms cancel), and the open-chain dispersion sum identity Σ_{k=1}^{⌊N/2⌋} ε(k)² = N−1 gives ‖M‖²_F = 4c²·(N−1)·4^N = 4·‖H‖²_F·2^N (with ‖H‖²_F = c²·(N−1)·2^N), exactly the F49 chain norm quoted in the Zero-Is-The-Mirror section above. M is γ-independent (Master Lemma); the norm scales with the coupling c, not γ.

- **Lebensader principle**: F80 is the third-layer manifestation of the Π·L·Π⁻¹ + L + 2σ·I = 0 funnel:
  - State layer: `cockpit_panel` (16 Paulis → 3-class trichotomy)
  - Real-space single-body operator layer: F78 (any (c, P) with P ∈ {Y, Z} → same M_l)
  - Momentum-space chain 2-body operator layer: F80 (4 Π²-odd Pauli pairs → same M-spectrum via Bloch sign-walk)

- **Π²-even XY-summed cousin at N=4** (operator-space Majorana lens witness): the same JW Bloch dispersion ε(k) = 2J·cos(πk/(N+1)) = {±φ, ±1/φ} governs the Π²-even chain H = (J/2)·Σ(XX+YY). Its axis-mode (n_XY = 2) Im(λ) clusters decompose into integer combinations of {φ, 1/φ, 1, √5}, and site-reflection R sorts them into R-parity-protected groups (±√5 R-even-only, ±1 R-odd-only, 18 silent R-even). The 18 silent modes are the operator-space Majorana self-conjugate sector at N=4. See [`experiments/MAJORANA_AXIS_MODES.md`](../../experiments/MAJORANA_AXIS_MODES.md).

---

## Numerical Verification Summary

| Test | Path | N range | Result |
|------|------|---------|--------|
| Bit-exact cluster prediction | `pi2_odd_universality_data_sweep.py` | 3-6 | All match |
| Bit-exact cluster prediction (N=7 full SVD) | `n7_bloch_signwalk_verification.txt` | 7 | All 4 clusters at predicted values, mult 4096 each |
| Independent eigsh check (N=7) | `n7_eigsh_check.txt` | 7 | Top eigenvalues of M·M† match SV² predictions |
| Universality across (X,Y)/(X,Z)/(Y,X)/(Z,X) | sweep | 3-6 | 100% (all 4 give bit-identical clusters) |
| Pytest lock | `test_F80_bloch_signwalk_chain_pi2_odd` | 4, 5 | Passes |

---

## Status

The open-chain two-body sign-walk locations are proven for the scoped
bilinears and verified at N=3–7. Their equal multiplicity is table-specific;
the collision-count formula above is the general statement. The direct
per-bond Π-action has wider finite checks, but topology-specific and k-body
cluster-value formulas are not established here.
