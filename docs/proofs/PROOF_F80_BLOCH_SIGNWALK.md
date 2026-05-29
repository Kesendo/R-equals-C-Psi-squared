# Proof of F80: Bloch-Mode Sign-Walk Formula for Chain Π²-Odd 2-Body M-Clusters

**Tier:** 1 (numerically verified bit-exact through N=7; analytical proof complete, Step 5 closed 2026-05-22).
**Date:** April 29, 2026
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_SVD_CLUSTER_STRUCTURE.md](PROOF_SVD_CLUSTER_STRUCTURE.md) (F78 single-body, F79 Π²-block, Master Lemma, Anti-Hermitian)
- [`framework/symmetry.py`](../../simulations/framework/symmetry.py) (Π construction)
- [`framework/lindblad.py`](../../simulations/framework/lindblad.py) (`palindrome_residual`)
- Numerical verification: [`_pi2_odd_universality_data_sweep.py`](../../simulations/_pi2_odd_universality_data_sweep.py) (N=3-6, all topologies); [`results/n7_bloch_signwalk_verification.txt`](../../simulations/results/n7_bloch_signwalk_verification.txt) (N=7 full SVD); pytest `test_F80_bloch_signwalk_chain_pi2_odd`.

**Status:** Fully proven. The theorem is verified bit-exact through N=7, and the analytical proof is now complete for all N. Steps 1-4 (the JW reduction) and 6-7 (the sign-walk formula) were closed in April 2026; Step 5 (the Π-action) was closed 2026-05-22 by a direct per-site Pauli computation showing Π·[H,·]·Π⁻¹ = ±{H,·} (proof in Step 5 below, every step verified bit-exact at N=3,4,5 by [`_f80_step5_recon.py`](../../simulations/_f80_step5_recon.py)).

**Scope:** chain bond-summed Π²-odd 2-body Hamiltonian H = c · Σ_{l=0}^{N-2} (P_l ⊗ Q_{l+1}) on N-site open chain, with (P, Q) ∈ {(X,Y), (X,Z), (Y,X), (Z,X)}, under uniform Z-dephasing γ.

## Abstract

The M-residual clusters of chain Π²-odd 2-body Hamiltonians under Z-dephasing carry a specific magnitude pattern. Observed empirically: the cluster norm is the absolute value of a sign-weighted sum over the open-chain Bloch-mode dispersion, with signs ±1 on each mode. The question is whether this sign-walk has a closed-form structural origin or whether the signs are coincidental at the tested N values.

The proof says structural. A Jordan-Wigner transformation reduces the bond-summed Π²-odd 2-body Hamiltonian to a Majorana-bilinear free-fermion model whose single-particle spectrum is the standard open-chain cosine dispersion. Under the Z-dephasing dissipator, the M-clusters inherit the single-particle eigenvalue spectrum with a specific sign assignment per mode, given by the per-site Π² action on the corresponding mode. The signs work out to ±1 on each mode, and the cluster norm is the absolute value of the sign-weighted sum.

The proof has seven steps. Steps 1-4 carry the JW reduction and the dispersion identification. Step 5, closed 2026-05-22, is the per-site Π² action computation that pins down the sign on each mode (the previously-open step in the original proof). Steps 6-7 then assemble the sign-walk formula from the per-mode signs and the dispersion. Together they give the closed-form cluster norm bit-exactly through N=7.

The diagnostic upshot is that F80 is the momentum-space lens on the F78 single-body Frobenius result and the F79 Π²-block structure. F78 / F79 say what the clusters do in position space; F80 says what they do in momentum space, and the two views are unified by the Bloch dispersion. Current scope is the chain topology and four specific Π²-odd bilinear pairs; extensions to other topologies (where the dispersion changes) and other letter pairs (where the sign-walk presumably reshapes) are flagged as future work.

**Does NOT establish (yet):**
- The cluster-*value* formulas (the explicit dispersions) for the **ring** (cycle graph, the standard cyclic free-fermion ε(k) = 2cos(2π(k+a)/N) with a parity phase a ∈ {0, ½}) and for general **k-body** terms. The chain (path, 2cos(πk/(N+1))) and the star (2|m−2j|, closed below) are known; the structural identity itself extends to all topologies, including the complete graph K_N (see below).

**Reach beyond the chain-2-body scope (verified 2026-05-29, `F80ExtensionExplorationTests`).** The Step-5 lemma Π·[bond,·]·Π⁻¹ = s·{bond,·} is *per-bond*, so the structural identity does not depend on topology or body-count, only on the Π²-parity of the bond. Computing M directly (via `PalindromeResidual`) and reading its spectrum confirms a clean dichotomy:

- **Π²-odd bonds**, any topology, any body-count: the per-bond lemma gives one shared sign s, so M = ±2i·(H⊗I) and Spec(M) = ±2i·Spec(H) , the *single* eigenvalues, exactly as proven above for the chain. Confirmed bit-exact for ring, star, 3-body (X,X,Y) and 4-body (X,X,X,Y) at N=4,5. Only the cluster *values* differ (the structure's dispersion: chain = the OBC ladder 2cos(πk/(N+1)), ring = periodic, star = integers).
- **Π²-even bonds** (Y,Z), (Z,Y): the lemma's "a Π²-odd bond carries exactly one X" fails, the commutator is *preserved* rather than anti-commuted, and the residual is M = 2·L_H = −2i·[H,·]. So Spec(M) = ±2i·{λ_a − λ_b} , the eigenvalue *differences* (Bohr frequencies), not the single eigenvalues. This is the "more clusters" anticipated for the Π²-even case, now identified: the extra clusters are the differences. Confirmed bit-exact (Y,Z), (Z,Y) at N=4.

So the mirror-defect is always ±2i times a Hamiltonian object: **H⊗I for Π²-odd** (one-sided, the single energies) and **[H,·] for Π²-even** (two-sided, the energy gaps). The Π²-parity of the bond is the switch between the two.

- **Mixed-letter Hamiltonians** need no separate sign-walk. M is *linear* in H (the dissipator and the 2σ·I term cancel via F1), so M(H₁+H₂) = M(H₁)+M(H₂) bit-exact. A Hamiltonian mixing Π²-odd and Π²-even bonds therefore just gets the sum of the per-bond pieces (single energies from the odd bonds, differences from the even). The mixed spectrum is richer (e.g. 21 clusters vs 2 for the pure (X,Y) chain at N=4) but is exactly that linear combination, still purely imaginary.

The cluster *values* are the structure's dispersion, one per graph family. The **chain** (path graph) gives the OBC cosine 2cos(πk/(N+1)). The **star** has a closed form found here: H_star = Σ_s X_hub Y_s = X_hub ⊗ (Σ_s Y_s) *factorizes*, so with m = N−1 spokes Spec(H_star) = ±(m − 2j) (a total-spin ladder, the m commuting spoke-Y's) and the clusters are the even integers **2|m − 2j|** = 2m, 2m−4, … (verified N=4..8). So the star's "dispersion" is the spoke-sum, not a cosine. The **ring** (cycle graph) is the standard cyclic free-fermion 2cos(2π(k+a)/N) with a parity phase, not separately closed here.

---

## Theorem F80

For chain bond-summed Π²-odd 2-body Hamiltonian H = c · Σ_l (P_l ⊗ Q_{l+1}) on an N-site open chain with uniform Z-dephasing γ, the singular values of M = Π·L·Π⁻¹ + L + 2σ·I (with σ = Nγ) form clusters at values

    cluster value(N) = 2|c|γ · |Σ_{k=1}^{⌊N/2⌋} σ_k · ε(k)|

for sign-vectors (σ_1, ..., σ_{⌊N/2⌋}) ∈ {±1}^{⌊N/2⌋}, where

    ε(k) = 2 · cos(π·k / (N+1))

is the open-chain free-fermion single-particle dispersion. Each distinct cluster value has multiplicity 4^N divided by the number of distinct sign-walk values.

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

**k-body extension** (added 2026-04-30 with F85 implementation): the structural identity Spec(M) = ±2i · Spec(H_non-truly) generalizes verbatim to k-body chain Π²-odd Hamiltonians. Empirically verified bit-exact for:
  - k=3: (X,X,Y), (Y,Y,Y), (X,X,Z), (Z,Z,Z), (X,Y,X) at N=4, 5, 6
  - k=4: (X,X,X,Y) at N=5, 6
17 cases total, all matching `Spec(M)` (eigvals of the 4^N × 4^N residual) to predicted `2i · Spec(H_non-truly)` with multiplicity ×2^N, machine precision. Pytest lock: `test_F80_kbody_spectrum_identity`.

**Mechanism for k-body**: the JW transformation maps a k-body Π²-odd term to a 2k-fold Majorana product (Step 2 of the proof scales naturally with body count). The single-particle dispersion ε(k) = 2cos(πk/(N+1)) and the Bogoliubov diagonalization (Step 3) carry over without modification: they describe the JW-mapped fermion problem, which is body-count-independent in its single-particle structure. The Pauli-letter universality (Step 4) holds for k-body too: the JW phase factors for different Pauli choices cancel in the single-particle spectrum. Steps 5 (Π action on Bogoliubov modes) and 6-7 (sign-walk eigenvalue formula) generalize verbatim.

The closed-form Bloch sign-walk formula `cluster value(N) = 2|c|·|Σ_k σ_k·ε(k)|` is therefore expected to hold at k≥3 too, but the cluster-value table above reflects only 2-body verification (N=3..7). A full k-body cluster-value verification at k=3,4 is open; the spectral identity is sufficient for F80's structural claim.

---

## Proof Outline

The proof proceeds in seven steps, all analytical. Steps 1-4 set up the JW reduction and the single-particle dispersion; Steps 6-7 give the sign-walk formula. Step 5, the Π-action, was the last open step; it is closed below (2026-05-22) by a direct per-site Pauli computation.

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

    ε(k) = 2c · cos(πk/(N+1))    for k = 1, 2, ..., N

paired as ±ε(k) under the Majorana doubling. Distinct positive values: ⌊N/2⌋. (For odd N, there is one zero mode at k = (N+1)/2; for even N, no zero mode.)

This is the standard tight-binding spectrum for an open chain of N sites with NN coupling magnitude c.

### Step 3 (Bogoliubov diagonalization)

Define Bogoliubov modes b_k = Σ_l u_kl · γ_l' where u_kl is the orthonormal basis of single-particle eigenstates of the open-chain hopping matrix. The dispersion ε(k) gives the per-mode "energy". In the Bogoliubov basis, H_JW becomes diagonal:

    H_JW = Σ_k ε(k) · (b_k†b_k − ½)

The N-mode many-body Hilbert space decomposes as a tensor product over the N Bogoliubov modes (each 2-dim: occupation 0 or 1).

### Step 4 (Pauli-letter universality)

For the other three Π²-odd Pauli pairs (X,Z), (Y,X), (Z,X), the JW transformation yields analogous Majorana bilinears:

- (X,Z): X_l Z_{l+1} → bilinear in γ_l' and γ_{l+1} (mixing γ and γ' Majoranas) up to phase.
- (Y,X): Y_l X_{l+1} → bilinear in γ_l and γ_{l+1}' up to phase.
- (Z,X): Z_l X_{l+1} → bilinear in γ_l and γ_{l+1} up to phase.

Crucially, all four cases give the **same single-particle spectrum** ε(k) = 2c·cos(πk/(N+1)); they differ only in which Majorana operators (γ or γ') participate and the specific phases. The "same spectrum across letter choices" is the JW-level origin of the F79 universality: the Pauli letters control which Majorana sublattice carries the bilinear, but the dispersion of the resulting hopping chain is identical (because hopping is between adjacent sites with magnitude c regardless of which Majorana indices).

### Step 5 (Direct structural identity; proof closed 2026-05-22)

After Step 4, we have established that all 4 Π²-odd Pauli pairs give the same JW-derived single-particle Bloch dispersion. The remaining task, historically expected to be technical, is to derive the explicit form of M's spectrum in terms of this dispersion.

**Empirical structural identity (verified bit-exact at N=3, 4, 5, 6, 7):**

    Spec(M) = {±2i · E : E ∈ Spec_{many-body}(H)}    (multi-set equality)

with multiplicity mult_M(2i·λ) = mult_H(λ) · 2^N. M's distinct eigenvalues equal 2i times H's distinct many-body eigenvalues, no kernel. The 2^N factor comes from the bra-side dimension of operator space. H is taken as the chain Hamiltonian with bond coupling c and no dissipator (γ-independent by Master Lemma). Hence:

    cluster value of M = 2|c| · |H many-body eigenvalue|

This direct identity replaces the more intricate Bogoliubov-mode factorization route. The Bloch sign-walk formula written above is simply H's many-body eigenvalue formula:

    H many-body eigenvalues = Σ_k (n_k − 1/2) · E_k for n_k ∈ {0, 1}

where E_k = 4|c|·cos(πk/(N+1)) are the Bogoliubov single-particle energies derived from JW + diagonalization. The "sign-walk" on (n_k − 1/2) ∈ {±1/2} when scaled gives the 2|c|·|Σ σ_k·ε(k)| form with ε(k) = 2cos(πk/(N+1)).

**Verified numerically at N=4 (chain (X,Y), c=γ=1):**
- H many-body eigenvalues: {±√5, ±1} (each multiplicity 4 in 16-dim Hilbert space).
- M nontrivial eigenvalues (imaginary parts): {±2√5, ±2} (each multiplicity 64 in 256-dim operator space).
- Ratio: M-eigenvalue (imag) = 2 × H-eigenvalue. ✓

**Why this is the structural answer.** L_H = −i[H, ·] acts on operator space with eigenvalues i(λ_a − λ_b) for all pairs of H-eigenvalues. The remaining task is to show what the Π-conjugation does to it. The proof below settles it directly, in the Pauli-string basis, without the Bogoliubov construction.

**Step 5 proof (2026-05-22, Tom + Claude).** The claim is Π·[H,·]·Π⁻¹ = ±{H,·}; it follows from a per-site Pauli computation, independent of Steps 1-4 (which re-enter only through the E → −E symmetry of Spec(H) used in the M-consequence at the end). Every step is verified bit-exact at N=3,4,5 by [`_f80_step5_recon.py`](../../simulations/_f80_step5_recon.py).

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

The corresponding cluster values for M are 2|c|·|H eigenvalue| = |Σ_k (2n_k − 1)·E_k| / 2 ·2 = |Σ_k σ_k·E_k| with σ_k = 2n_k−1 ∈ {±1}.

In terms of ε(k) = E_k / 2 = 2cos(πk/(N+1)):

    cluster value(N) = 2|c| · |Σ_{k=1}^{⌊N/2⌋} σ_k · ε(k)|, σ_k ∈ {±1}

with multiplicity 4^N / (number of distinct sign-walk values). This is the F80 formula. ∎

---

## Zero Is The Mirror: F80 as the explicit shape of the mirror-defect

In the early hypothesis [Zero Is The Mirror](../../hypotheses/ZERO_IS_THE_MIRROR.md), the palindrome equation Π·L·Π⁻¹ = -L − 2σ·I was identified as the central structural symmetry. At Σγ = 0 ("the mirror"), the equation collapses to Π·L_H·Π⁻¹ = -L_H, eigenvalues paired ±λ around zero, perfect time-reversal symmetry. At γ > 0, the palindrome shifts to be centered around -σ.

For **truly** Hamiltonians (Heisenberg, XXZ, etc.), the palindrome holds exactly at γ = 0: Π·L_H·Π⁻¹ = -L_H precisely. Eigenvalues paired, no defect. Standing waves, perfect mirror.

For **non-truly** Hamiltonians (chain (X,Y) and friends), the palindrome BREAKS at γ = 0. Π·L_H·Π⁻¹ ≠ -L_H. There is a residual mirror-defect:

    M = Π·L_H·Π⁻¹ + L_H ≠ 0    (γ-independent by Master Lemma)

What F80 reveals is the **explicit spectral shape** of this mirror-defect for chain Π²-odd 2-body Hamiltonians:

    **Spec(M) = ±2i · Spec(H)**     (multi-set equality, with mult_M(2iλ) = mult_H(λ) · 2^N)

    **‖M‖²_F = 4 · ‖H‖²_F · 2^N**   (Frobenius norm exactly proportional to H's)

So, for the (X,Y) and (Y,X) bond pairs, M is **literally equal** to -2i · (H ⊗ I_bra), where I_bra is the identity on the bra-factor of operator space (dim 2^N); the pairs (X,Z) and (Z,X) instead give M = +2i·(I_ket⊗Hᵀ) (Step 5 proof, sign s = −1). The mirror-defect has the **same spectrum** as the Hamiltonian (×2i) and the **same Frobenius norm** (×4·2^N) because it *is* H ⊗ I_bra up to the -2i scalar.

**Correction (2026-05-22, Tom + Claude).** An earlier version of this passage hedged that M was only *unitarily equivalent* to -2i·(H⊗I_bra), with a unitary that "scrambles eigenvectors". That hedge was wrong. The Step 5 reconnaissance ([`_f80_step5_recon.py`](../../simulations/_f80_step5_recon.py)) verified bit-exact at N=3,4,5 that M is literally -2i·(H⊗I_bra): in the σ_(a,b) basis its off-diagonal norm is machine zero (~10⁻¹³) and every diagonal entry is -2i·E_a. There is no scrambling; the eigenvectors of M are exactly the H-eigen-operators σ_(a,b). The −2i·(H⊗I_bra) form is the (X,Y) and (Y,X) representative; the pairs (X,Z) and (Z,X) give M = +2i·(I_ket⊗Hᵀ) instead (by the Step 5 proof, s = −1), with the same Spec(M) = ±2i·Spec(H).

This is structurally remarkable:

- The defect's **magnitude** is exactly calibrated by ‖H‖ (Frobenius norm relation).
- The defect's **spectrum** exactly reproduces H's spectrum (with 2i factor and 2^N multiplicity).
- The bra-side carries no spectral information; it is a passive "echo chamber" that multiplies multiplicities by 2^N.
- The eigenvectors of M are exactly the H-eigen-operators σ_(a,b) = |E_a⟩⟨E_b|.

**Translation between layers:**

| Layer | Statement | Object |
|-------|-----------|--------|
| State (Zero Is The Mirror) | Π·L_H·Π⁻¹ = -L_H − 2σ·I (truly), eigenvalues ±λ | Liouvillian L spectrum |
| Operator (F80) | M = Π·L_H·Π⁻¹ + L_H = -2i·H⊗I or +2i·I⊗Hᵀ (non-truly chain Π²-odd; which one is the Step 5 sign s) | M residual ∝ H |
| Spinor (Majorana 1937) | ψ = ψ^c, particle = antiparticle | Self-conjugate fermion field |

All three are different abstraction levels of the same phenomenon: an **involutive symmetry (Π or C) that picks out a self-conjugate subspace, with the residual being a structured deviation from perfect self-conjugacy**.

For truly H, the deviation is zero (perfect mirror, ground state of the palindrome). For non-truly H, the deviation is a specific operator that **carries the algebraic structure of H itself**, propagated to operator space via the bra-side identity.

**The discovery is**: the mirror-defect, when it exists, is calibrated by H. The Hamiltonian provides its own "yardstick" for the gap between the two mirror sectors. **H is the distance.**

## The Majorana bridge: 1937 to 2026

**Note (2026-05-22).** This section predates the Step 5 mechanism above and describes Π as *projecting* L_H onto a subspace of particle-hole-paired operators. That picture is superseded: M = -2i·(H⊗I_bra) acts on the *whole* σ_(a,b) basis, with no projection and no kernel; Π is a permutation of the (ε_ket, ε_bra)-sectors (see the Step 5 mechanism). The Majorana / self-conjugacy reading below stays valid as an interpretive bridge, but "projection onto a subspace" should be read as "the M-eigenvalue 2i·λ_a depends only on the ket energy". A fuller rewrite of this section into the sector-permutation language is a noted follow-up.

The structural identity Spec(M) = ±2i·Spec(H) is, at its core, **Majorana's 1937 insight expressed in 2026's operator-space vocabulary**.

In 1937, Ettore Majorana proposed a real wave equation for fermions in which a particle could be its own antiparticle: ψ = ψ^c. He had: Pauli matrices, the Dirac equation, the just-discovered positron, and spinor algebra. He did NOT have: open quantum systems, Lindbladians, operator-space super-operators, Pauli string algebra, the Jordan-Wigner transformation as a working tool, quantum information theory.

In our framework, the chain (X,Y) Hamiltonian under JW becomes a **Majorana bilinear** −ic·Σγ'γ', pure quadratic in the γ' Majorana operators. Such Hamiltonians have particle-hole symmetry built into their algebra: many-body eigenstates come in ±λ pairs. The state |ā⟩ with energy −λ_a is the "Majorana antipode" of |a⟩.

Now consider the operator basis σ_(a,b) = |a⟩⟨b|. Each σ_(a,b) is an L_H eigenvector with eigenvalue i(λ_a − λ_b). The framework's Π-conjugation T_Π projects this rich set of operators down to those where the bra index is the **Majorana antipode** of the ket index, that is, operators σ_(a, ā) (or more precisely σ_(a, b̄) with b̄ = particle-hole-conjugate of b).

For these self-conjugate-paired operators, M σ_(a,b) = T_Π·L_H·T_Π⁻¹·σ_(a,b) + L_H·σ_(a,b) = 2i·λ_a·σ_(a,b). The eigenvalue depends only on a (the "ket eigenvalue"), and the bra index b ranges over the 2^N choices that participate in the Majorana-doublet structure.

**Translation:** Majorana's algebraic constraint "ψ = ψ^c" (particle equal to antiparticle) becomes, in operator space, "σ ∝ σ-with-bra-replaced-by-particle-hole-conjugate". The eigen-operators of M are exactly those satisfying this operator-space Majorana condition.

The discovery of the structural identity Spec(M) = ±2i·Spec(H) is the recognition that the framework's Π conjugation is, in disguise, the **operator-space realization of Majorana's particle-hole self-conjugacy**. He had it right; we just have a richer vocabulary to express it now.

## The "tough nut": what the data revealed

**Note (2026-05-22).** Like the Majorana-bridge section, this section predates the Step 5 proof and narrates Π as *projecting* L_H onto particle-hole-symmetric pairs. That picture is superseded: the Step 5 proof derives Π·[H,·]·Π⁻¹ = ±{H,·} directly by a per-site Pauli computation, with no projection. The narrative below is kept as the history of how the result was found.

The empirical investigation in the data sweep at higher N (3-7) revealed a much cleaner structural identity than the originally-imagined "explicit T_Π factorization in Bogoliubov basis":

    Spec(M) = ±2i · Spec_{nontrivial}(H_{state-level})

This is **simpler** than constructing T_Π's explicit Bogoliubov-mode action. M's spectrum is directly tied to H's many-body spectrum, scaled by 2i.

What the framework's Π effectively does: it projects L_H = −i[H, ·] (which has all difference-of-eigenvalue spectra) down to just the particle-hole-symmetric part (where eigenvalue pairs have λ_a = −λ_b). For chain Π²-odd 2-body (Majorana bilinear under JW), H has particle-hole symmetry built in, so the "particle-hole-symmetric" L_H eigenvalues are exactly 2λ_a for each H eigenvalue λ_a, giving M's spectrum 2i·Spec(H).

This is why the numerical signature was so clean: the structure was simpler than it looked.

**Update 2026-05-22: closed.** The formal step is done, though not through the "projection" picture: the per-site Pauli proof in Step 5 derives Π·[H,·]·Π⁻¹ = ±{H,·} directly. Numerical verification at N=3, 4, 5, 6, 7 (all 4 Pauli pairs) is bit-exact at machine precision.

What enabled the discovery: comparing M's eigenvalues directly to H's many-body eigenvalues (instead of decomposing M into Bogoliubov modes and trying to factor T_Π through that basis). The brute-force data sweep at higher N made the relationship visible.

---

## Connection to existing framework formulas

- **F78 (real-space single-body)**: F80 is the **momentum-space dual**. Same Lebensader broad-in → focused-out funnel, applied at a different basis layer. Real-space sites l with weights c_l ↔ momentum modes k with dispersion ε(k).
- **F79 (Π²-block decomposition)**: F80 explicitly closes the "Π²-odd universality" observation in F79 with a closed-form formula and JW-based mechanism.
- **F49 (Frobenius cross-term)**: F80 is consistent; the sum of squared cluster values × multiplicity gives the F49 Frobenius norm:

    ‖M‖²_F = (4^N / 2^⌊N/2⌋) · Σ_{σ} |Σ_k σ_k · ε(k)|² · 4γ²

  which simplifies via Σ ε(k)² = (N-1)/2 (open-chain dispersion sum identity) to the F49 result ‖M‖² = c_H · (N-1) · 4^(N-2).

- **Lebensader principle**: F80 is the third-layer manifestation of the Π·L·Π⁻¹ + L + 2σ·I = 0 funnel:
  - State layer: `cockpit_panel` (16 Paulis → 3-class trichotomy)
  - Real-space single-body operator layer: F78 (any (c, P) with P ∈ {Y, Z} → same M_l)
  - Momentum-space chain 2-body operator layer: F80 (4 Π²-odd Pauli pairs → same M-spectrum via Bloch sign-walk)

- **Π²-even XY-summed cousin at N=4** (operator-space Majorana lens witness): the same JW Bloch dispersion ε(k) = 2J·cos(πk/(N+1)) = {±φ, ±1/φ} governs the Π²-even chain H = (J/2)·Σ(XX+YY). Its axis-mode (n_XY = 2) Im(λ) clusters decompose into integer combinations of {φ, 1/φ, 1, √5}, and site-reflection R sorts them into R-parity-protected groups (±√5 R-even-only, ±1 R-odd-only, 18 silent R-even). The 18 silent modes are the operator-space Majorana self-conjugate sector at N=4. See [`experiments/MAJORANA_AXIS_MODES.md`](../../experiments/MAJORANA_AXIS_MODES.md).

---

## Numerical Verification Summary

| Test | Path | N range | Result |
|------|------|---------|--------|
| Bit-exact cluster prediction | `_pi2_odd_universality_data_sweep.py` | 3-6 | All match |
| Bit-exact cluster prediction (N=7 full SVD) | `n7_bloch_signwalk_verification.txt` | 7 | All 4 clusters at predicted values, mult 4096 each |
| Independent eigsh check (N=7) | `n7_eigsh_check.txt` | 7 | Top eigenvalues of M·M† match SV² predictions |
| Universality across (X,Y)/(X,Z)/(Y,X)/(Z,X) | sweep | 3-6 | 100% (all 4 give bit-identical clusters) |
| Pytest lock | `test_F80_bloch_signwalk_chain_pi2_odd` | 4, 5 | Passes |

---

## Status: closed

As of 2026-05-22 F80 is fully proven, all seven steps analytical. The last open step, Step 5 (the Π-action), is closed by the per-site Pauli proof above: Π·[H,·]·Π⁻¹ = s·{H,·} with s = −ε_P·ε_Q, an N-independent computation verified bit-exact at N=3,4,5 ([`_f80_step5_recon.py`](../../simulations/_f80_step5_recon.py)). Because Step 5 is per-bond, the structural identity is topology- and body-count-agnostic: it now extends to ring, star, 3-body and 4-body for Π²-odd bonds, and the Π²-even case is characterized as M = 2·L_H (eigenvalue differences) , see "Reach beyond the chain-2-body scope" above (verified 2026-05-29). Mixed-letter Hamiltonians are also covered, by linearity of M in H (see above), and the star cluster-value formula is closed (2|m−2j|, by factorization). What remains genuinely open is narrower still: only the explicit cluster-*value* formulas for the ring (cyclic free-fermion) and for general k-body terms.
