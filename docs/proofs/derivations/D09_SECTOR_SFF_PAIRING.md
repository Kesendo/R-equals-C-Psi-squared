# D09: Sector SFF Pairing

**Derives:** K_freq(w, t) = K_freq(N−w, t) for all t
**From:** Formula 1 (palindrome equation) + Π weight complementarity
**Status:** PROVEN (corollary of mirror symmetry proof)

---

## Statement

The spectral form factor (SFF) restricted to XY-weight sector w is
identical to the SFF of sector N−w:

    K_freq(w, t) = K_freq(N−w, t)    for all t ≥ 0

where K_freq(w, t) = |Σ_k exp(i·ω_k·t)|² summed over eigenvalues
λ_k = −d_k + i·ω_k in sector w.

## Proof

**Step 1. Π maps sector w to sector N−w.**

The conjugation operator Π maps each Pauli string σ with XY-weight
w_XY(σ) = w to a string with XY-weight N − w. This is proven in
Step 1 of the [Mirror Symmetry Proof](../MIRROR_SYMMETRY_PROOF.md):
Π swaps {I,Z} ↔ {X,Y} at each site, so every non-XY factor becomes
XY and vice versa.

**Step 2. Π induces a bijection on eigenvalues within paired sectors.**

Since Π L Π⁻¹ = −L − 2Σγ·I (the palindrome identity), applying Π to
an eigenvector v_k with eigenvalue λ_k = −d_k + i·ω_k in sector w
yields an eigenvector Π·v_k in sector N−w with eigenvalue
λ_k' = −(2Σγ − d_k) − i·ω_k.

This map is a bijection: every eigenvalue in sector w has a unique
partner in sector N−w, and vice versa.

**Step 3. Paired sectors have identical frequency sets (up to sign).**

The eigenvalue map gives:
- d_k' = 2Σγ − d_k (palindromic pairing of decay rates)
- ω_k' = −ω_k (frequency negation)

Therefore the frequency multiset of sector N−w is {−ω_k : λ_k ∈ sector w}.

**Step 4. SFF is invariant under frequency negation.**

The SFF is defined as:

    K_freq(w, t) = |Σ_k exp(i·ω_k·t)|²

Replacing ω_k → −ω_k:

    K_freq(N−w, t) = |Σ_k exp(−i·ω_k·t)|²
                   = |conj(Σ_k exp(i·ω_k·t))|²
                   = |Σ_k exp(i·ω_k·t)|²
                   = K_freq(w, t)

since |z|² = |z̄|² for any complex number z.    ∎

## Special case: XOR sector

The XOR sector (w = N) maps to w = 0 (all I/Z strings). Both sectors
consist entirely of real eigenvalues (ω = 0), so K_freq = (count)²
for both. At w = N, all eigenvalues are degenerate at rate 2Nγ
(Formula 23), giving K = (N+1)² × δ(t). The w = 0 sector has the
N+1 stationary + near-stationary modes with the palindromic partner
rates 2Σγ − 2Nγ = 2(Σγ − Nγ).

## Verification

Numerically verified for N=3−5 (sector-by-sector SFF comparison,
relative error < 10⁻¹²). See [`simulations/spectral_form_factor.py`](../../../simulations/spectral_form_factor.py).
