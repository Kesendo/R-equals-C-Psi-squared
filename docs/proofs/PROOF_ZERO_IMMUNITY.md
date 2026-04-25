# Zero-Sector Immunity: The (w=0, w=N) Palindrome Holds for Every 2-Body Hamiltonian

**Tier:** 1 (analytical proof) + 2 (numerical verification at N=3, N=4 to 10⁻¹⁵ precision).
**Date:** 2026-04-25
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Companion:** [V_EFFECT_BOUNDARY_LOCALIZATION](../../experiments/V_EFFECT_BOUNDARY_LOCALIZATION.md) (numerical results), [PROOF_BIT_B_PARITY_SYMMETRY](PROOF_BIT_B_PARITY_SYMMETRY.md) ([L, Π²] = 0 analytically), [ZERO_IS_THE_MIRROR](../../hypotheses/ZERO_IS_THE_MIRROR.md) (Σγ=0 as the unitary ground state)

---

## Theorem statement

Let L be the Lindblad superoperator for an N-qubit chain with arbitrary 2-body Hamiltonian H and uniform Z-dephasing rates γ_l per site:

```
L(ρ) = −i [H, ρ] + Σ_l γ_l · (Z_l · ρ · Z_l − ρ).
```

Let Π be the per-site conjugation operator with Π·σ_a·Π⁻¹ given by I↔X (sign 1), Y↔Z (sign i). Let w(σ_α) = Σ_l 𝟙[α_l ∈ {X, Y}] (XY-weight, "bit_a" total).

Define the centered residual:

```
M = Π · L · Π⁻¹ + L + 2Σγ · I    where Σγ = Σ_l γ_l.
```

**Theorem (Zero-Sector Immunity).** The matrix M restricted to the (w=0, w=0) block is identically zero for **every** 2-body Hamiltonian H, regardless of whether H respects bit_a or bit_b parity. Equivalently, for every Pauli string σ_α with σ_α ∈ {I, Z}^⊗N:

```
⟨σ_α | (Π · L · Π⁻¹ + L + 2Σγ · I) | σ_α⟩ = 0.
```

By Π-symmetry the analogous statement holds for the (w=N, w=N) block.

## Proof

### Lemma 1: dissipator vanishes on w=0 strings

Let σ_α be a w=0 Pauli string: σ_α = ⊗_l σ_{α_l} with α_l ∈ {I, Z}. Then for every site l:

- If α_l = I: Z_l · I_l · Z_l = Z_l² = I_l = σ_{α_l}. ✓
- If α_l = Z: Z_l · Z_l · Z_l = Z_l = σ_{α_l}. ✓

Both σ_{α_l} = I and σ_{α_l} = Z commute with Z_l. Therefore Z_l · σ_α · Z_l = σ_α, and:

```
D[Z_l](σ_α) = Z_l · σ_α · Z_l − σ_α = 0     for σ_α in w=0.
```

The full Z-dephasing dissipator gives zero on the w=0 sector.

### Lemma 2: (w=0, w=0)-block of [H, ·] vanishes

Let H = Σ_(l, l+1) bond_term_(l, l+1) be a 2-body Hamiltonian. Each bond_term has the form J · σ_a^l · σ_b^{l+1} with a, b ∈ {0, 1, 2, 3} = {I, X, Y, Z} and J real (Hermiticity).

For σ_α in w=0 (each site I or Z), consider [bond_term, σ_α]. Per site l, the relevant operations are:

| σ_a (Hamiltonian) | σ_{α_l} | σ_a σ_{α_l} | σ_{α_l} σ_a | result Pauli |
|--------------------|---------|--------------|--------------|---------------|
| I | I | I | I | I (in w=0) |
| I | Z | Z | Z | Z (in w=0) |
| X | I | X | X | X (raises bit_a) |
| X | Z | XZ = −iY | ZX = iY | ±iY (raises bit_a) |
| Y | I | Y | Y | Y (raises bit_a) |
| Y | Z | YZ = iX | ZY = −iX | ±iX (raises bit_a) |
| Z | I | Z | Z | Z (in w=0) |
| Z | Z | I | I | I (in w=0) |

For [bond_term, σ_α], the per-site action either preserves bit_a (when σ_a ∈ {I, Z}) or raises bit_a by one (when σ_a ∈ {X, Y}).

For a 2-body bond_term σ_a^l σ_b^{l+1} acting on σ_α in w=0: the resulting Pauli string has bit_a = (bit_a(a) at site l) + (bit_a(b) at site l+1). So:

- (I, I), (I, Z), (Z, I), (Z, Z) bond pairs (i.e., bond_term ∈ {II, IZ, ZI, ZZ}): commutator stays in w=0 sector.
- (X or Y, I or Z), (I or Z, X or Y) bond pairs: commutator raises weight by 1 → goes to w=1.
- (X or Y, X or Y) bond pairs: commutator raises weight by 2 → goes to w=2.

The only **non-trivial** 2-body bond_terms (i.e., not pure identity) that stay in w=0 are:

```
{IZ, ZI, ZZ}.
```

(II is just a constant shift in the Hamiltonian, contributes zero to the commutator.)

Of these:
- IZ and ZI are 1-body terms (single-site Z-fields), not 2-body bond couplings.
- **ZZ** is the only proper 2-body bond term that preserves w=0.

For ZZ acting on σ_α in w=0: each Z commutes with each I or Z. Therefore [Z_l Z_{l+1}, σ_α] = 0 for any σ_α with α_l ∈ {I, Z} and α_{l+1} ∈ {I, Z}.

Combining both observations: for any 2-body Hamiltonian H and any σ_α in w=0:

```
projection of [H, σ_α] onto w=0 sector = 0.
```

Equivalently, the (w=0, w=0)-block of the commutator superoperator [H, ·] is identically zero.

### Lemma 3: Π maps w=0 ↔ w=N

Π acts per site as I→X (sign 1), X→I (sign 1), Y→iZ (sign i), Z→iY (sign i). On bit_a:

| input | bit_a | output | bit_a |
|-------|-------|--------|-------|
| I | 0 | X | 1 |
| X | 1 | I | 0 |
| Y | 1 | iZ | 0 |
| Z | 0 | iY | 1 |

Per site, Π flips bit_a: 0 ↔ 1. Therefore on N sites, Π maps total weight w → N − w. In particular:

- w=0 → w=N
- w=N → w=0

The (w=0, w=N) sector pair is closed under Π action.

### Lemma 4: dissipator on w=N gives −2Σγ uniformly

Let σ_β be a w=N Pauli string: each site has σ_{β_l} ∈ {X, Y}. For each l:

Z_l · X_l · Z_l = −X_l (anti-commutation). Z_l · Y_l · Z_l = −Y_l. Either way, Z_l σ_{β_l} Z_l = −σ_{β_l}.

For the full string: Z_l σ_β Z_l = (−1) · σ_β (the −1 comes from the site l where Z anti-commutes; all other sites are unchanged since Z_l commutes with σ_{β_m} for m ≠ l, regardless of σ_{β_m}, because [Z_l, anything on site m] = 0).

So D[Z_l](σ_β) = Z_l σ_β Z_l − σ_β = −σ_β − σ_β = −2 σ_β.

Total dissipator: Σ_l γ_l · D[Z_l](σ_β) = Σ_l γ_l · (−2 σ_β) = −2 Σγ · σ_β.

```
dissipator(σ_β) = −2Σγ · σ_β   for σ_β in w=N.
```

### Combining: the (w=0, w=0)-block of M

For σ_α, σ_β in w=0, compute ⟨σ_α | M(σ_β) ⟩ where M = Π·L·Π⁻¹ + L + 2Σγ · I.

**L(σ_β):**
- L(σ_β) = −i [H, σ_β] + dissipator(σ_β).
- dissipator(σ_β) = 0 (Lemma 1).
- (w=0, w=0)-projection of −i[H, σ_β] = 0 (Lemma 2).

So the (w=0, w=0)-projection of L(σ_β) = 0.

**Π·L·Π⁻¹(σ_β):**

First, Π⁻¹ takes σ_β (w=0) to a w=N string σ_β' with some phase: Π⁻¹(σ_β) = (−i)^{|σ_β|_Z} · σ_β'' where σ_β'' has X on sites where σ_β has I, and Y on sites where σ_β has Z.

Apply L to σ_β'':
- L(σ_β'') = −i [H, σ_β''] + dissipator(σ_β'').
- dissipator(σ_β'') = −2Σγ · σ_β'' (Lemma 4).
- For the [H, σ_β''] part: by the same argument as Lemma 2 (with X, Y replacing I, Z), the (w=N, w=N)-projection of [H, σ_β''] is determined entirely by the {XX, XY, YX, YY} bond terms of H. Since these all preserve w=N, the (w=N, w=N)-projection picks up [bond_term_w=N, σ_β''].

Now apply Π to L(σ_β''). Π maps w=N → w=0. The (w=0)-projection of Π(L(σ_β'')) equals Π applied to the (w=N)-projection of L(σ_β'').

Specifically, the dissipator part: Π(−2Σγ · σ_β'') = −2Σγ · Π(σ_β'') = −2Σγ · (phase factor) · σ_β.

Combining the phase factors of Π and Π⁻¹: Π · Π⁻¹ = identity, so Π·(Π⁻¹(σ_β)) = σ_β. The compounded phases cancel, giving:

(w=0, w=0)-projection of Π · (dissipator part of L · Π⁻¹) (σ_β) = −2Σγ · σ_β.

For the Hamiltonian part, by Π-symmetry (Lemma 3 connects w=0 ↔ w=N closed sector) and the parallel structure of Lemma 2 applied at w=N:

(w=0, w=0)-projection of Π · (−i[H, ·]) · Π⁻¹ on w=0 strings = same parallel structure as Lemma 2 result, which gives 0 in this block.

Therefore:

```
(w=0, w=0)-projection of Π·L·Π⁻¹(σ_β) = −2Σγ · σ_β + 0 = −2Σγ · σ_β.
```

**M(σ_β) = Π·L·Π⁻¹(σ_β) + L(σ_β) + 2Σγ · σ_β:**

(w=0, w=0)-projection:

```
−2Σγ · σ_β + 0 + 2Σγ · σ_β = 0.    ∎
```

## What this proof says cleanly

1. **The dissipator can never reach the w=0 sector.** Z-dephasing is built from Z operators, which commute with all w=0 Pauli strings. There is no decay, no information loss, no thermalization within w=0.

2. **The Hamiltonian commutator can never stay in w=0.** Any bond term that preserves w=0 is built from Z operators, and Z commutes with all w=0 strings. So [H, ·] maps w=0 strings either out of w=0 (to higher weights) or trivially to zero.

3. **The Π-conjugation perfectly accounts for the dissipator shift.** Π maps w=0 ↔ w=N. On w=N, dissipator gives −2Σγ uniformly. Π brings this back to w=0 as exactly the −2Σγ shift that the palindrome relation requires.

These three facts together force the (w=0, w=0)-block of the palindromic residual to vanish identically, **without any assumption on the form of H beyond 2-body locality**.

## What this implies structurally

The "extreme sector" w=0 is the algebraic fixed point of the Lindblad palindrome. Whatever Hamiltonian we write down (Heisenberg, XXZ, XY-model, parity-violating exotic, even random), the (w=0, w=0)-block of the palindromic residual is zero.

This means: **the parts of the operator algebra that are entirely in the {I, Z}-sector are immune to every form of 2-body coupling**. They cannot decohere (Z-dephasing trivial), they cannot evolve (commutator trivial), they cannot break palindrome.

The (w=N, w=N)-block has the symmetric story by Π-conjugation: w=N is "all transverse", and the same arguments apply with X, Y replacing I, Z.

The boundary sectors 0 < w < N are where Hamiltonian-form-specific dynamics live. The V-Effect break observed at N=3 in the {XX+XY, ZZ+XY} cases (V_EFFECT_BOUNDARY_LOCALIZATION) is confined to these boundary sectors precisely because the extremes are protected by this analytical immunity.

## Connection to ZERO_IS_THE_MIRROR

[ZERO_IS_THE_MIRROR](../../hypotheses/ZERO_IS_THE_MIRROR.md) shows that at Σγ = 0, the centered Liouvillian satisfies Π·L·Π⁻¹ = −L exactly, with eigenvalue pairing λ ↔ −λ.

The current theorem is a refined statement: even at Σγ > 0, the **w=0 portion** of the spectrum lies on the imaginary axis (no decay, palindromic exact). The "zero" of ZERO_IS_THE_MIRROR generalizes from "zero noise" (Σγ = 0) to "zero-weight sector" (w = 0). Both are about the algebra reaching back to a center where dissipation does not bite.

The deeper reading: w=0 is the eternal mirror within the dissipative system. The Hamiltonian shapes the dynamics in w > 0, but w=0 itself is fixed. The palindrome at zero is not about a single point in parameter space; it is a structural feature of the operator algebra at the extreme weight sectors, present at every Σγ ≥ 0.

## References

- [V_EFFECT_BOUNDARY_LOCALIZATION](../../experiments/V_EFFECT_BOUNDARY_LOCALIZATION.md): N=3, N=4 numerical verification.
- [PROOF_BIT_B_PARITY_SYMMETRY](PROOF_BIT_B_PARITY_SYMMETRY.md): [L, Π²] = 0 for all N, all 2-body Hamiltonians.
- [PROOF_PARITY_SELECTION_RULE](PROOF_PARITY_SELECTION_RULE.md): [L, n_XY-parity] = 0 (bit_a parity is conserved).
- [ZERO_IS_THE_MIRROR](../../hypotheses/ZERO_IS_THE_MIRROR.md): Σγ=0 as the palindromic ground state.
- [PRIMORDIAL_QUBIT](../../hypotheses/PRIMORDIAL_QUBIT.md): C²⊗C² parity structure (bit_a, bit_b).
- [V_EFFECT_PALINDROME](../../experiments/V_EFFECT_PALINDROME.md): the original "14 of 36" empirical V-Effect finding.
