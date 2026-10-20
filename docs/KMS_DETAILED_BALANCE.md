# KMS and Detailed Balance: Locating Π Among Liouvillian Symmetries

<!-- Keywords: Pi operator not KMS detailed balance, shifted anti-similarity
Liouvillian, quantum detailed balance Alicki 1977, Roberts hidden time-reversal
symmetry PRX Quantum 2021, tenfold Lindbladian classification Sa Prosen 2023,
Buca-Prosen weak strong symmetry gap, pure Z-dephasing algebra,
finite-temperature thermal-jump comparison, 2:2 split Pi operator not palindrome,
algebraic palindrome shift
2Sgamma, R=CPsi2 KMS detailed balance -->

**Status:** KMS/QDB comparison complete; full irreducible-sector SRP class OPEN; a finite-temperature element of Π's shape CONFIRMED for XXZ bonds with axial DM terms, without fields or Z dephasing, where one ratio holds across the hopping: the cooling-to-heating ratio, or with X or Y dephasing the one-site steady state's, whose square root is then the weight; no site-local element exists where that ratio differs across a hopping bond
**Last refreshed:** 2026-10-20
**Authors:** Thomas Wicht, Claude (Anthropic, Cowork Research)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md), [Π as Time Reversal](../experiments/PI_AS_TIME_REVERSAL.md)

---

## What this document is about

Physics has well-established notions of "time-reversal symmetry" for
open quantum systems: the KMS condition (from statistical mechanics)
and quantum detailed balance (from thermodynamics). This document asks
whether our conjugation operator Π is one of those known symmetries.
The Π anti-similarity is not itself a KMS or quantum-detailed-balance
condition. That does not imply that the same generator fails every QDB
definition: Alicki/standard QDB permits a separate Hamiltonian derivation,
and at the faithful invariant state `I/d` the Heisenberg-plus-Z-dephasing
generator satisfies that weaker decomposition. In the Sá-Ribeiro-Prosen (SRP)
framework, however, the centered generator and irreducible symmetry sectors
are the correct objects: Π gives a P-type anticommutation there after resolving
Π². The final SRP class of our fully reduced sectors has not been computed.
The pure-Z channel algebra alone does not identify a thermal bath or a temperature.
A thermal bath does give a mirror of Π's shape for XXZ bonds with axial DM terms
without fields or dephasing, built from its weight where one cooling-to-heating
ratio holds across the hopping, and with X or Y dephasing from the weight of the
one-site steady state where that state's ratio does (Question 4).

---

## Abstract

The Π relation is **not** a quantum-detailed-balance (QDB) or KMS condition and
is not a standard Buca-Prosen symmetry. It is a shifted anti-similarity:
Π·L·Π⁻¹ = −L − 2Sγ·I. Strict GNS/KMS symmetry of the full generator relates
L to its weighted adjoint and implies a real spectrum; Π instead relates L to
its negative and gives palindromic complex pairs. Alicki/standard QDB permits a
separate Hamiltonian derivation, so only the dissipative part is self-adjoint
and that condition does not force the full generator to have a real spectrum.
At β=0, strict full-generator KMS symmetry reduces to Hilbert-Schmidt
self-adjointness and fails whenever H≠0, while the weaker Alicki decomposition
can still hold. The closest existing framework is Roberts-Lingenfelter-Clerk
hidden time-reversal symmetry (PRX Quantum 2021), which also produces
time-reversal-like structure in systems violating standard DB. The SRP
classification is formulated for a trace-shifted generator, so the constant
2Sγ is absorbed by `L_c = L + SγI` rather than creating a classification gap.
Globally Π has order four because Π² is a commuting unitary symmetry; within
each Π²-parity sector a phase-normalized restriction is an involutive P generator.
Assigning BDI, CI, or another SRP class still requires the complete symmetry
algebra after all strong/unitary symmetries are resolved. The known Π construction uses the 2:2
per-site Pauli split of pure Z dephasing. That is an algebraic property of the
jump channel, not a bath-temperature assignment. The spectral palindrome for
thermal jump operators is a separate question; for a Hamiltonian that commutes with X⊗N and Y⊗N it
survives at the centre −Σ(γ↓+γ↑)/2 ([F137](ANALYTICAL_FORMULAS.md)), and where one
cooling-to-heating ratio holds across the hopping, for XXZ bonds with axial DM terms,
no fields and no dephasing, an exact site-local element of Π's shape carries it,
built from the thermal weight; with X or Y dephasing the one-site steady state takes
the bath's place (Question 4). Here `2Sγ` is the algebraic
shift of the pure-Z palindrome, not an entropy-production rate.

---

## Executive Summary

The Π anti-similarity is **not** a standard quantum detailed balance (QDB)
condition. Its closest relative in the literature is a recently identified class
of **hidden time-reversal symmetries** in open quantum systems
(Roberts-Lingenfelter-Clerk, PRX Quantum 2021), which also gives time-reversal-like
structure where standard detailed balance fails; the formal connection is not
established (Question 6). The comparison
to KMS at β=0 is formal, not a bath assignment: a generator that is strictly
KMS-symmetric as a whole at equilibrium `I/d` is self-adjoint (real spectrum),
while Π gives
palindromic spectral pairing (complex spectrum with μ → -μ). These are different
symmetries with different consequences. For SRP classification one first centers
L, resolves the commuting Π² sectors, and phase-normalizes Π in each sector. The
remaining sectorwise class is OPEN until the full symmetry algebra is computed.
At a finite temperature, a site-local element of Π's shape exists for XXZ bonds
with axial DM terms, without fields or dephasing, where one cooling-to-heating
ratio holds across the hopping; it is built from the thermal
weight of detailed balance, and with X or Y dephasing from the one-site steady
state's (Question 4).

---

## Question 1: Is Π a Known Symmetry Type?

**Classification: CONFIRMED algebra / OPEN class and prior-art equivalence**

### The Buca-Prosen framework (2012)

Buca and Prosen classified Liouvillian symmetries as:

- **Weak symmetry:** [U, L] = 0 (unitary U commutes with the Liouvillian)
- **Strong symmetry:** U·L_k = e^{iθ}·L_k·U for each jump operator L_k individually

Our Π satisfies neither. It does not commute with L. Instead:

    Π · L · Π⁻¹ = -L - 2Sγ · I

This is an **anti-similarity** with a constant shift. If we define L_c = L + Sγ·I
(the centered Liouvillian), then Π·L_c·Π⁻¹ = -L_c, which is a pure anti-commutation.
In the Buca-Prosen language, Π is a "weak anti-symmetry" of L_c, but this category
does not exist in their 2012 classification. They only considered [U, L] = 0.

### The Albert-Jiang framework (2014)

Albert and Jiang extended the Buca-Prosen classification to include conserved
quantities and the structure of the steady-state manifold. Their framework
characterizes the infinite-time behavior of Lindblad evolution, including
degenerate steady states. However, their classification focuses on operators
that commute with L (symmetries), not operators that anti-commute (anti-symmetries).
They do mention anti-unitary symmetries in passing, but do not develop the theory
for the anti-commuting case.

### The tenfold Lindbladian classification (Sá et al., PRX 2023)

The most comprehensive classification is the **38-fold symmetry classification**
of many-body Lindbladians by Sá, Ribeiro, and Prosen (Phys. Rev. X 13, 031019, 2023).
This classifies Lindbladians by their behavior under:

- Two flavors of time-reversal: T₊ and C₊
- Two flavors of particle-hole: T₋ and C₋
- Chiral/sublattice symmetry: P
- Pseudo-Hermiticity: Q₊
- Anti-pseudo-Hermiticity: Q₋

SRP applies negative symmetries to the trace-shifted generator. For this model
that object is `L_c = L + SγI`, on which Π anticommutes exactly. Moreover
Π² = U_X commutes with L_c. In a U_X-parity sector p_x = ±1,
`P_{p_x} = sqrt(p_x) Π|_{p_x}` satisfies `P_{p_x}² = I` and is therefore an
ordinary sectorwise P generator. The full class requires the remaining SRP
symmetries and their commutation signs in each irreducible sector.

Their local-dephasing examples land in BDI or CI depending on the resolved parity
sector and Hamiltonian. Those labels are a source-level precedent, not labels we
can copy onto our Hamiltonian without the corresponding sector calculation.

### Assessment

The established statement is algebraic: Π is an invertible superoperator that
implements a shifted anti-similarity of L and a P-type anticommutation of L_c.
Its global order four is resolved sectorwise through Π² = U_X. No global AIII,
BDI, or CI label is assigned here; the irreducible-sector classification remains
open.

**Literature and scope comparison.** We built this from the dephasing algebra with
no literature input: the −2Sγ shift fell out of the 2:2 per-site Pauli split, and the
shifted-sublattice form below was recognized afterward, on 2026-06-08. The pairing about a nonzero center,
λ ↔ −λ − 2Sγ (eigenvalue pairs ±λ + ia), has the shifted-sublattice-symmetry
form catalogued by Kawasaki, Mochizuki, and Obuse (Phys. Rev. B 106, 035408,
2022): ordinary chiral/sublattice symmetry plus a constant decay-rate shift.
Their construction is noninteracting at the quadratic/Majorana level and uses
`S = I₂ ⊗ τ_z`. The object analyzed here is the shifted anti-similarity for the
stated interacting/k-body model together with a bounded 6 → 4 → 2
per-site-product-versus-entangled locality result. Prior-art coverage and
equivalence for that interacting, locality-resolved scope remain OPEN. This is
an object-and-scope comparison, not an ownership, priority, or novelty claim.

---

## Question 2: Is Π Related to Quantum Detailed Balance?

**Classification: DISTINCT CONDITIONS for Π; the thermal element M_β of Question 4 is built from the detailed-balance weight**

### Detailed-balance conventions that must not be conflated

For a faithful invariant state `ρ`, let `L̃` denote the adjoint in the chosen
`ρ`-weighted inner product. A **strict GNS/KMS-symmetry** convention imposes

    L̃ = L.

That strict condition makes the full generator self-adjoint in that inner
product. Alicki/standard QDB is weaker: in a common Heisenberg-picture
convention one decomposes

    L = L₀ + i[K, ·],    [K, ρ] = 0,    L̃₀ = L₀.

Thus the reversible Hamiltonian derivation is allowed and only the dissipative
part `L₀` is self-adjoint. Conditions involving an antiunitary time reversal or
a modular transform are related detailed-balance definitions, not automatically
equivalent to `L̃=L` without their additional hypotheses.

Here `L̃` is the adjoint for the GNS product tr(ρ a* b) or the KMS product
tr(ρ^½ a* ρ^½ b); both reduce to the Hilbert-Schmidt product at ρ = I/d, so `L̃`
becomes the Hilbert-Schmidt adjoint `L†` there (Question 3). `ρ_ss` below is the
same ρ, the steady state that serves as the weight.

**What QDB says physically:** the forward and backward transition rates between
any two states are related by the Boltzmann factor. The Liouvillian "looks the
same" under time reversal weighted by the thermal state.

**Spectral boundary:** strict full-generator GNS/KMS symmetry implies real
eigenvalues. Alicki/standard QDB does not force the full generator to have a
real spectrum, because `i[K,·]` can supply oscillatory imaginary parts.

### Our Π condition

    Π · L · Π⁻¹ = -L - 2Sγ · I

**What Π says:** L is similar to its own negative (up to a shift). Every
eigenvalue λ has a partner -(λ + 2Sγ).

**What Π implies spectrally:** eigenvalues come in palindromic pairs with
BOTH real parts (decay rates) and imaginary parts (frequencies) paired.
The spectrum is complex, not real.

### Why they are different

The key structural differences:

| Property | Standard QDB | Our Π condition |
|---|---|---|
| Relates L to | strict form: L̃ = L; Alicki form: dissipative L₀ is self-adjoint and a Hamiltonian derivation is separate | -L - c (negative + shift) |
| Sign | Same in the strict symmetric part | Opposite (minus sign) |
| Constant shift | None | 2Sγ |
| Spectral consequence | strict full-generator form: real; Alicki form: full spectrum may be complex | Palindromic complex pairs |
| Involves adjoint? | Yes (the weighted adjoint L̃) | No (L itself) |
| Steady state role | Central (ρ_ss is the weight) | None; the shift is algebraic |

Detailed balance constrains a weighted adjoint, either for the full generator
in the strict symmetry convention or for its dissipative part after a commuting
Hamiltonian derivation is split off. Our Π is a statement about L being similar
to its own *negative*. These are distinct constraints.

The table compares Π. At a finite temperature the weight does enter a mirror: the
element of Π's shape for a thermal bath, M_β of Question 4, is F137's S_XY followed
by the transpose of the state and the detailed-balance weight of
[PROOF_ANTILINEAR_TRIANGLE](proofs/PROOF_ANTILINEAR_TRIANGLE.md) §7, the square root
of the stationary state on both sides.

**A Lindbladian can satisfy both constraints, one, or neither.** Pure dephasing
without a Hamiltonian is Hilbert-Schmidt self-adjoint and has a real palindromic
spectrum; that algebraic fact does not choose a thermal bath. Heisenberg plus
all-site Z dephasing has the faithful invariant state `I/d`; its dissipator is
Hilbert-Schmidt self-adjoint and its Hamiltonian is an allowed derivation, so it
satisfies Alicki/standard QDB with respect to `I/d` while failing strict
full-generator self-adjointness when `H≠0`. Independently, it satisfies Π and
can have a complex palindromic spectrum.

### The Alhambra-Woods connection (2017)

Alhambra and Woods ("Dynamical maps, quantum detailed balance, and the Petz
recovery map," 2017) prove a hypothesis-scoped connection to recovery maps.
For a QDB dissipative semigroup with no unitary part and a full-rank invariant
state, the dynamical map is equal to its Petz recovery map. With an additional
commuting unitary part, the Petz map reverses the unitary sign while retaining
the same dissipative evolution. This is not a blanket equivalence between every
QDB convention and exact reversal of the full channel. Our Π does not have
either recovery interpretation: it pairs modes rather than reversing the
unitary evolution or reproducing the dissipative map.

---

## Question 3: Does KMS at β=0 Reduce to Our Π?

**Classification: NO (different symmetries; the pure-Z channel has no β)**

### The β=0 KMS reference point

For a generator that is strictly KMS-symmetric as a whole and whose faithful
equilibrium state is
`ρ_ss = I/d`, the modular operator is

    Δ = ρ_ss ⊗ ρ_ss⁻¹ = (I/d) ⊗ (dI) = I.

Its strict symmetry condition therefore becomes

    L† = L     (self-adjointness in Hilbert-Schmidt norm).

This is a statement about strict full-generator KMS symmetry with a specified
equilibrium state, not about every convention called quantum detailed balance
and not a consequence of writing down a pure-Z jump operator.

### The pure-Z model does not select β = 0

The pure-Z dissipator is unital and fixes every density matrix diagonal in the
Z basis. It does not alone select `I/d`, a unique steady state, or a Gibbs β.
With a Hamiltonian, its stationary manifold is determined by the full
generator: the selected finite unbiased-TFI runs reach `I/d`, whereas the
all-site-dephased Heisenberg branch retains the F4 sector structure. Neither
case makes the pure-Z channel a KMS bath at β = 0.

For the Hamiltonian-plus-dephasing generator,

- `L_H(ρ) = -i[H, ρ]` is anti-self-adjoint: `L_H† = -L_H`;
- `L_D` is self-adjoint in the Pauli basis with real eigenvalues; and
- `L = L_H + L_D` therefore has `L† = -L_H + L_D ≠ L` whenever `H ≠ 0`.

So the strict β=0 KMS-symmetry condition is not Π and is not implied by the
pure-Z model. The separate Alicki/standard-QDB decomposition can nevertheless
hold at `I/d`, because it retains the anti-self-adjoint Hamiltonian derivation.

### What the model does satisfy

Instead it satisfies

    Π · L · Π⁻¹ = -L - 2Sγ · I.

This is a different symmetry. It does not require self-adjointness, permits
complex eigenvalues, and pairs rates algebraically. The shift `2Sγ` contains no
Gibbs state or bath temperature.

---

## Question 4: Does a Finite-Temperature Generalization Exist?

**Classification: CONFIRMED for XXZ bonds with axial DM terms, without fields or Z dephasing, where the one-site steady-state ratio is the same on every set of sites joined by hopping (an exact site-local element of Π's shape, built from the thermal weight, and with X or Y dephasing from the steady state's); no site-local element where that ratio differs across a hopping bond; OPEN beyond it**

**Swept for the element.** [ANALYTICAL_FORMULAS](ANALYTICAL_FORMULAS.md): F137 holds the
thermal palindrome and the element S_XY of the shape −L† − c, and its one-ratio carrier
holds the weighted transpose as an equality of spectra; F119 and F84 hold the transpose
law and the β = 0 reading; F108 Part 2 holds the X-dephasing Π, Π_5b(X), per site
I → Z, Z → −I, X → −iY, Y → iX, and
[PROOF_F108_KLEIN_V4_EQUIVALENCE](proofs/PROOF_F108_KLEIN_V4_EQUIVALENCE.md) its other
orientation Π_5b(X)⁻¹; `PiOperator` types a further phase variant, and
[MIRROR_SYMMETRY_PROOF](proofs/MIRROR_SYMMETRY_PROOF.md) records that either sign of the
X ↔ Y phase works.
[PROOF_ANTILINEAR_TRIANGLE](proofs/PROOF_ANTILINEAR_TRIANGLE.md) §7 holds the weight
itself, L(H)ᵀ = Γ_β⁻¹·L(Hᵀ)·Γ_β, and says it composes with S_XY into the palindrome,
again for spectra; [PROOF_PI_FACTORS_AS_R_TIMES_D](proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md)
holds the transpose inside Π_Z. [THE_ONE_SQUARE](THE_ONE_SQUARE.md) §9,
MIRROR_SYMMETRY_PROOF and [THERMAL_BREAKING](../experiments/THERMAL_BREAKING.md) hold
S_XY. OQ-168 in [the open-questions index](../review/OPEN_QUESTIONS_INDEX.md) holds the
question of the thermal operator. MirrorWorld's `AntilinearTriangle` pins the weighted
transpose; its `Sextant` searches the one-site closure, which holds no weight, and finds
S_XY times the copies and nothing of Π's shape in a thermal bath at ratio ⅓. The OpenArcs
registry, `fw.Confirmations`, the glossary, [CAUGHT_ERRORS](CAUGHT_ERRORS.md), the Core
claims and the Diagnostics witnesses: nothing on a weighted mirror. No store held a
single operator of the shape −L − c at finite temperature. Checked for adjacency: the
weight against §7's (the same object); the letter map below at equal rates against
F108 Part 2's X-dephasing Π (per site −i·τ is Π_5b(X)⁻¹, the other orientation);
τ against Π_Z (τ(ρ) = Y⊗N·Π_Z(ρ), the shift by Y⊗N after Π_Z, two reflections of the
Z-dephasing dissipator, so τ keeps that dissipator, which is why Z dephasing breaks
the element).

### The rates

At finite temperature the bath has, on each site, the jump operators

    J↓ = √(γ(n̄+1)) σ⁻    (emission, rate γ↓ = γ(n̄+1))
    J↑ = √(γn̄) σ⁺          (absorption, rate γ↑ = γn̄)

where n̄ = 1/(exp(ℏω/kT) − 1) is the Bose-Einstein occupation number and
σ⁻ lowers the qubit.

The emission/absorption ratio is γ↓/γ↑ = (n̄+1)/n̄ = exp(ℏω/kT).
At T→∞: ratio → 1, giving symmetric excitation and relaxation; this is not
pure Z dephasing.
At T→0: ratio → ∞, pure emission (amplitude damping).

Our Π works because the per-site Pauli rates have a 2:2 split
({I,Z} immune, {X,Y} decaying at 2γ). At finite T, the per-site
Pauli rates for the thermal bath are:

- I: rate 0 (always immune)
- X: rate γ(2n̄+1)/2 (nonzero, but less than Z)
- Y: rate γ(2n̄+1)/2 (same as X)
- Z: rate γ(2n̄+1) (nonzero; Z is no longer immune!)

The rates are [0, r/2, r/2, r] where r = γ(2n̄+1) = γ↓ + γ↑. This is NOT the
1:3 isotropic split of depolarizing noise. The rate-pairing condition
from DEPOLARIZING_PALINDROME.md requires the 4 rates to partition
into 2 equal-sum pairs. Here: (0, r) and (r/2, r/2) both sum to r.
So rate-pairing IS possible in principle.

### Where the pairing has to sit

Rate-pairing alone is not sufficient: the element must also
anti-commute with [H, ·] (Step 2 of the proof). Our known Π
maps {I,Z} ↔ {X,Y}, which is exactly the swap needed for Z-dephasing
where {I,Z} are immune and {X,Y} decay.

For the thermal bath rates [0, r/2, r/2, r], the rate-pairing requires
I ↔ Z (rates 0 and r) and X ↔ Y (rates r/2 and r/2). This pairing does not
swap diagonal and off-diagonal sectors; it permutes within each. Written
as a multiplication on both sides, ρ ↦ Y⊗N·ρ·X⊗N, it sends I to −iZ, Z to iI,
X to Y and Y to X on each site. A two-sided map ρ ↦ A·ρ·B with
A·H·A⁻¹ = B⁻¹·H·B = H′ carries [H, ·] to [H′, ·]; the strings Y⊗N and X⊗N
leave real XXZ bonds unchanged, so this one commutes with [H, ·] there and
cannot give −L − c. Through the transpose the same pairing does: the transpose
reverses products, [H, ρ]ᵀ = −[Hᵀ, ρᵀ], and

    τ(ρ) = Y⊗N · ρᵀ · X⊗N

(I to −iZ, Z to iI, X to Y, Y to −X on each site) anticommutes with [H, ·]
exactly when X⊗N·H̄·X⊗N = H (H̄ the entry-wise complex conjugate) and H commutes
with Z⊗N, since then
X⊗N·Hᵀ·X⊗N = Y⊗N·Hᵀ·Y⊗N = H. XXZ bonds with axial DM terms on any graph,
real hopping included and fields excluded, satisfy it. Per site, −i·τ is the
inverse of F108 Part 2's X-dephasing Π, the other orientation that
PROOF_F108_KLEIN_V4_EQUIVALENCE records; the inverse of a reflection reflects too.
For depolarizing noise even
rate-pairing is impossible (1:3 split); for the thermal bath it is possible, and
the transpose carries it past the Hamiltonian.

### The element

Notation: N̂ counts the excitations; c = Σ_l (γ↓_l + γ↑_l) is the shift;
S_XY is F137's element ρ ↦ X⊗N·ρ·Y⊗N, which gives −L† − c for an H commuting
with X⊗N and Y⊗N; D[J](ρ) = J·ρ·J† − ½{J†J, ρ}; σ⁻ = |0⟩⟨1| lowers the excited
state |1⟩; Γ_β is ρ ↦ π^½·ρ·π^½ for the weight π defined below.

At equal cooling and heating, γ↑_l = γ↓_l on every site, τ itself reflects L,
τ·L·τ⁻¹ = −L − c, and dephasing along X or Y beside the damping is reflected
with it, as under the X-dephasing Π. At a ratio other than 1 the jumps carry an
arrow (the generator applied to the identity is (γ↓ − γ↑)·Z on each site), and
τ needs the thermal weight. Let π ∝ e^{−βN̂} with e^{−β} = γ↑/γ↓, the stationary
product state of the damping (β in units of the bare qubit splitting, as in
PROOF_ANTILINEAR_TRIANGLE §7), and let the ratio γ↑_l/γ↓_l be the same on every
set of sites joined by hopping, the rates themselves free; on each such set π
takes that set's β, so π is a product over the sites. Then

    M_β(ρ) = π^½ · Y⊗N · ρᵀ · X⊗N · π^½,      M_β · L · M_β⁻¹ = −L − c.

M_β is linear and site-local, a product of one-site maps, since the transpose, the
strings and the weight all factor over the sites. The proof has three lines.
(i) On one site, M₁(ρ) = w·Y·ρᵀ·X·w with w = diag(1, q) and q² = γ↑/γ↓ carries
γ↓·D[σ⁻] + γ↑·D[σ⁺] to its negative minus (γ↓ + γ↑), for symbolic rate and q.
(ii) τ anticommutes with [H, ·], as above, and Γ_β commutes with [H, ·]
because the weight commutes with H: hopping conserves N̂ and joins only sites of
one β, and ZZ bonds are diagonal. (iii) The dissipator is a sum over sites, each
acted on by its own factor. The two conditions in (ii) are different: τ needs
X⊗N·H̄·X⊗N = H with [H, Z⊗N] = 0, which XYZ bonds meet without conserving N̂;
the weight needs the excitation number conserved; at a ratio of 1 the weight is
the identity and only the first is used.

What the element is made of. At q ≠ 1 it goes beyond Pauli permutations (w mixes
I with Z), and it carries a Boltzmann factor, of the bare energy N̂ rather than of
H, in the symmetric split π^½·(·)·π^½ that KMS detailed balance uses (detailed
balance in the K-twisted (Agarwal) sense of PROOF_ANTILINEAR_TRIANGLE §7, not the
Alicki form of Question 2; F84 and MIRROR_SYMMETRY_PROOF use the words
for the β = 0 point γ↓ = γ↑ only). As maps on ρ it is S_XY followed by the transpose of the state and the weight,
M_β = (−1)^N·Γ_β∘(ρ ↦ ρᵀ)∘S_XY. The pieces were in the repository: S_XY in F137,
the weight in
[PROOF_ANTILINEAR_TRIANGLE](proofs/PROOF_ANTILINEAR_TRIANGLE.md#7-one-level-up-the-generators-transpose-and-the-two-ends-of-a-channel)
§7, where it enters through the transpose of the generator; and followed by ρ ↦ ρ†,
M_β is, up to the same sign, the antilinear K∘Γ_β∘S_XY, with K the entry-by-entry
conjugation. Read the other way, it is the X-dephasing Π carried to a finite
temperature by the thermal weight.

Exact over ℚ(i) in [`cube_moves_gate.py`](../simulations/cube_moves_gate.py) M11:
the one-site identity with symbolic rate and q; the anticommutation of τ on a real
triangle and on one with hoppings 1, 2 and 3i at Δ = 7⁄10; the full element on that
triangle at γ↑ = γ↓/4 with three different rates; two ratios, a hopping pair
and a lone site joined to it by a ZZ bond; the ratio-1 map against the inverse of
F108 Part 2's Π, per site. The controls that must fail: the element without
the weight; a Z field on one site; Z dephasing (M_β commutes with its dissipator);
X dephasing beside the damping at γ↑ = γ↓/4; three ratios with the weight matched
site by site on a connected triangle, a weight that does not commute with the
hopping; the two components joined by a hopping bond instead.

Where it stops. A temperature gradient along the hopping, a Z field and Z
dephasing break the element; X or Y dephasing at a ratio other than 1 breaks it
only through its weight (below). Pure cooling (γ↑ = 0) and pure heating (γ↓ = 0) make the weight
singular; under cooling alone F137 carries the palindrome by its route through
the effective Hamiltonian ([F137](ANALYTICAL_FORMULAS.md)). A broken element is not
always a broken spectrum: with real hopping a gradient still pairs, through S_XY
(F137's paragraph on loss and gain exchanged), while with complex hopping a
gradient can break the pairing itself, and there no element of any kind exists.
With the weight of the one-site steady state in place of the bath's,
w = diag(1, q) with q² = (γ↑ + γ⊥)/(γ↓ + γ⊥) for dephasing at γ⊥ along any axis
in the XY plane, the same map reflects L onto −L − c, c = Σ(γ↓ + γ↑ + 2γ⊥) with
the dephasing included, wherever that ratio is finite, nonzero and the same on
every set of sites joined by hopping. Where it differs across a hopping bond no
site-local element of this shape exists at all: applied to the identity, such
an element would have to produce a product of one-site operators Z_l + κ_l·I,
κ_l = (q_l² − 1)/(q_l² + 1), commuting with H, and a hopping bond between
unequal κ does not ([THE_CUBE_SITE_BY_SITE](../experiments/THE_CUBE_SITE_BY_SITE.md)). Where fields break M_β,
F137's chained carrier still gives an element of the shape −L − c: S_XY, the
thermal transpose, the conjugation and the turns W about Z with W·(−H̄)·W† = H
are each an identity between generators, and composed they are one site-local
linear map, ρ ↦ W·π^½·X⊗N·ρ·Y⊗N·π^½·W†, wherever that chain applies (no
dephasing, one ratio on every set of sites joined by hopping, and at Δ = 0 turns
exist exactly when e^{2iΦ} = (−1)^ℓ on every cycle; [F137](ANALYTICAL_FORMULAS.md)).
Whether an element of the shape −L − c exists wherever the spectrum pairs and
neither M_β nor that chain carries it, a product over sites or not, is open, and
so is whether the hidden time reversal of Roberts et al. in the doubled space is
M_β read through the Choi-Jamiolkowski map.

### Assessment

The known Π requires {I,Z} to be immune (2:2 split), whereas the thermal jump
model gives Z a nonzero decay rate, so Π itself does not carry over. The
pairing the thermal rates ask for, I ↔ Z and X ↔ Y, does carry over, written
through the transpose and weighted by the thermal state: M_β gives Π's shape
−L − c exactly wherever one ratio holds across the hopping. Without the weight,
the phased pairing S_XY gives −L† − c at any rates for an H that commutes with
X⊗N and Y⊗N ([THE_ONE_SQUARE](THE_ONE_SQUARE.md) §9).

**The spectrum.** The Liouvillian spectrum of XXZ under a thermal bath is
palindromic exactly, about the centre −Σ(γ↓+γ↑)/2, at every N through S_XY
(above; [F137](ANALYTICAL_FORMULAS.md)), and checked by rational characteristic
polynomial at the rate profiles of
[`thermal_palindrome_centre.py`](../simulations/thermal_palindrome_centre.py) at
N = 2 and 3. The rate that appears is r = γ↓ + γ↑ of the
rates above. The operator of the shape −L − c at one
ratio is the element M_β above. Where a palindrome is observed without a known
element, construct Π numerically from eigenvector
pairing (the fallback construction in NON_HEISENBERG_PALINDROME.md; note it
returns an entangled representative under spectral degeneracy even when a
product Π exists, so test candidate products directly, see
PI_OPERATOR_ENTANGLEMENT.md).

---

## Question 5: What Is 2Sγ?

**Classification: CONFIRMED (algebraic shift)**

`2Sγ = 2 Σᵢ γᵢ` is the shift in the pure-Z palindrome. For the dissipator
alone, it is the decay rate of a Pauli string carrying `X` or `Y` at every
site. With a Hamiltonian, it remains the center-setting algebraic quantity in

    Π · L · Π⁻¹ = -L - 2Sγ · I.

It is not an entropy-production rate, a free-energy difference, or evidence
for a thermal equilibrium. Those meanings require a specified bath and a
selected stationary state, neither of which is supplied by pure Z dephasing.

---

## Question 6: Who Else Has Found Similar Structures?

### Roberts, Lingenfelter, Clerk (PRX Quantum 2, 020336, 2021)

**"Hidden Time-Reversal Symmetry, Quantum Detailed Balance and Exact
Solutions of Driven-Dissipative Quantum Systems"**

This is the closest existing work. Key parallels:

| Feature | Roberts et al. | Our Π |
|---|---|---|
| Type of symmetry | Hidden TRS (antiunitary) | Anti-similarity (linear) |
| Where it acts | Doubled system (thermofield double) | Single Liouvillian |
| Consequence | Time-symmetric correlation functions | Palindromic spectrum |
| Relates to detailed balance | Yes (defines a generalized DB) | Parallel but distinct |
| Works for driven systems | Yes | Only for dephasing + Hamiltonian |

The table compares Π. The thermal element M_β of Question 4, followed by ρ ↦ ρ† the antilinear
K∘Γ_β∘S_XY up to sign, carries the bath's weight; whether it is their hidden time reversal
read in the single Liouville space is open.

Roberts et al. show that hidden TRS provides exact steady-state solutions.
Our Π provides exact spectral pairing. Both are "time-reversal-like"
symmetries in systems that violate standard detailed balance.

**Key difference:** Their hidden TRS requires the doubled Hilbert space
and thermofield double state (a purification trick that represents a thermal mixed state as a pure state in a doubled system). Our Π acts directly on the Liouvillian
superoperator space. These might be related through vectorization (the
Choi-Jamiolkowski isomorphism, which maps superoperators to operators in the
doubled Hilbert space), but the formal connection is not established.

Their follow-up paper (PRL 134, 130404, 2025) extends hidden TRS to
driven spin chains and finds new dissipative phase transitions, suggesting
this is a rich framework with ongoing development.

### Sá, Ribeiro, Prosen (Phys. Rev. X 13, 031019, 2023)

**"Symmetry Classification of Many-Body Lindbladians: Tenfold Way and Beyond"**

Their classification is the most comprehensive framework for Liouvillian
symmetries. They identify 38 symmetry classes (10 without conserved
quantities, more with). They explicitly build dephasing examples.

Their negative symmetries are defined on the trace-shifted generator. Here the
centered Liouvillian `L_c = L + SγI` satisfies `Π·L_c·Π⁻¹ = -L_c`, so the
constant shift is already absorbed by the SRP construction. Because Π² = U_X
commutes with L_c, one first fixes U_X parity p_x; the phase-normalized
restriction `sqrt(p_x) Π` then squares to identity and supplies a P generator.
The final class is OPEN until the remaining symmetries and their algebra are
computed in fully irreducible sectors. SRP's own dephasing examples resolve to
BDI or CI depending on parity and Hamiltonian; those labels are not transferred
to this model without that calculation.

### Medvedyeva, Essler, Prosen (PRL 117, 137202, 2016)

**"Exact Bethe Ansatz Spectrum of a Tight-Binding Chain with Dephasing Noise"**

Their η-pairing symmetry in the Hubbard mapping is a related 1D free-fermion
construction. The compared scopes are:

- Their result is restricted to tight-binding (free fermion) chains in 1D
- The present algebra covers interacting spins (XXZ with arbitrary δ) on any graph
- Their η-pairing comes from the Hubbard model's SU(2) symmetry
- Here Π is defined by the Pauli algebra's {I,Z} ↔ {X,Y} structure

There is no explicit discussion of detailed balance in their paper.
The η-pairing is used to identify steady states and leading decay modes,
not to establish spectral pairing.

### Haga et al. (2023) - Incoherentons

Their incoherenton framework grades Liouvillian eigenmodes by XY-weight
(= incoherenton number). This is the natural language for our palindrome.
They see the band structure but do not identify the palindromic pairing
within or between bands. No discussion of detailed balance or KMS.

### Chen, Kastoryano, Gilyén (Comm. Math. Phys. 2025)

**"Efficient Quantum Gibbs Samplers with KMS Detailed Balance Condition"**

They construct Lindbladians satisfying a strict KMS-detailed-balance symmetry
for preparing thermal states. The symmetric generators studied there have real
spectra. Relevant for understanding what that strict KMS-DB condition looks
like mathematically, but their framework is for Gibbs
sampling (driving to thermal states), not for understanding dephasing
dynamics. No palindromic structures.

---

## Summary Table

| Question | Answer | Classification |
|---|---|---|
| 1. Is Π a known symmetry type? | For the centered generator `L_c=L+SγI`, Π supplies a sectorwise P generator after resolving Π²; the complete irreducible-sector SRP algebra is not yet computed | **CONFIRMED algebra / OPEN class** |
| 2. Is Π related to quantum detailed balance? | The Π anti-similarity is a distinct constraint. Strict GNS/KMS symmetry uses a weighted adjoint; Alicki/standard QDB may instead split off a commuting Hamiltonian derivation | **DISTINCT CONDITIONS** |
| 3. Does KMS at β=0 reduce to Π? | No. Strict full-generator KMS symmetry at `I/d` gives `L†=L`; Alicki/standard QDB can retain a Hamiltonian derivation. Neither condition reduces to Π | **NO** |
| 4. Finite-T element of **Π's shape**? | Π itself does not carry over (the 2:2 split is pure Z dephasing's). The thermal split [0, r/2, r/2, r] pairs as (I,Z) and (X,Y); as a multiplication on both sides that pairing commutes with [H, ·] for real XXZ bonds, through the transpose it anticommutes, and weighted by the thermal state it is the site-local linear element M_β(ρ) = π^½·Y⊗N·ρᵀ·X⊗N·π^½ with M_β·L·M_β⁻¹ = −L − c, c = Σ(γ↓ + γ↑), at one cooling-to-heating ratio across the hopping, for XXZ bonds with axial DM terms on any graph without fields or dephasing | **CONFIRMED** at one ratio (exact, `cube_moves_gate.py` M11), and with X or Y dephasing where the one-site steady-state ratio is the same along the hopping, weighted by that steady state ([THE_CUBE_SITE_BY_SITE](../experiments/THE_CUBE_SITE_BY_SITE.md)); a ratio that differs across a hopping bond leaves no site-local element; a Z field and Z dephasing break M_β, while with fields F137's chained carrier composes into a site-local element of another shape wherever its turns exist; beyond it **OPEN** |
| 4b. Finite-T generalization of the **palindrome**? | Not obstructed. The pairing the body computes here is real: the spectrum stays palindromic, centred at −Σ(γ↓+γ↑)/2. See [F137](ANALYTICAL_FORMULAS.md) | **CONFIRMED** (H = 0 by the tensor sum of the one-site eigenvalues {0, −r/2, −r/2, −r}, as F137's T1 tensor sum with r in place of γ, and the bottom line below; Heisenberg H exact at the rational rate profiles of thermal_palindrome_centre.py at N=2, 3 and measured to N=5, derived at every N for an H that commutes with X⊗N and Y⊗N by ρ ↦ X⊗N·ρ·Y⊗N, which gives −L† − c, [THE_ONE_SQUARE](THE_ONE_SQUARE.md) §9) |
| 5. What is 2Sγ? | Algebraic shift of the pure-Z palindrome; not a thermodynamic quantity without a specified bath and stationary state | **CONFIRMED** |
| 6. Who else has similar structures? | Roberts et al. (hidden TRS), Sá-Ribeiro-Prosen (shifted-generator classification and BDI/CI dephasing examples), and MEP (η-pairing) provide related structures; no identity or novelty claim is made without the remaining sector calculation | **CONFIRMED relatives / OPEN relation** |

---

## The Bottom Line

The Π relation is not a quantum-detailed-balance or KMS condition: it is a
**shifted anti-similarity of the Liouvillian** for the pure-Z-dephasing algebra,
and the algebra does not
itself identify a bath temperature. In the Sá-Ribeiro-Prosen framework the
shift is absorbed by `L_c=L+SγI`; after resolving Π², a phase-normalized Π is
a P generator in each parity sector. This places the structure inside the
shifted-generator classification problem rather than outside it. The final
BDI, CI, or other irreducible-sector label is **OPEN** until the remaining
unitary, strong, and antiunitary symmetries and their commutation algebra are
computed.

The closest living relative is the Roberts-Lingenfelter-Clerk hidden TRS
framework, which also produces time-reversal-like symmetry in systems that
violate standard detailed balance. Establishing the formal connection
(possibly via the Choi-Jamiolkowski mapping between superoperator space
and doubled Hilbert space) is the most promising direction for connecting
Π to the broader literature.

The thermal jump model has no 2:2 per-site split, so Π itself does not carry
over: that split is a property of pure Z dephasing, a channel algebra and not an
intrinsically infinite-temperature phenomenon. The rate pairing the thermal
model asks for does carry over. Written through the transpose and weighted by
the thermal state it is a site-local linear element of Π's shape at one
cooling-to-heating ratio, M_β(ρ) = π^½·Y⊗N·ρᵀ·X⊗N·π^½, with M_β·L·M_β⁻¹ = −L − c
(Question 4): the finite-temperature mirror is built from the thermal weight of
detailed balance in its K-twisted (Agarwal) form.

**Nor is there an obstruction to the palindrome.** The per-site thermal rates
[0, r/2, r/2, r], r = γ↓ + γ↑, pair as (0, r) and (r/2, r/2); their negatives are the
one-site eigenvalues, since the arrow I → Z is the generator's only off-diagonal
entry and makes it triangular, so the tensor sum over sites gives a palindrome
centred at −Σ(γ↓ᵢ + γ↑ᵢ)/2, and that survives the Heisenberg Hamiltonian too, exact at the script's rational rate profiles at N = 2 and 3, measured through 5, and derived at every N for an H that
commutes with X⊗N and Y⊗N by the element ρ ↦ X⊗N·ρ·Y⊗N, which gives −L† − c. See [F137](ANALYTICAL_FORMULAS.md). That
element needs no common ratio; M_β needs one across the hopping, and with real
hopping a gradient breaks M_β while the spectrum still pairs.

---

## Key References

1. **Alicki (1976)**: "On the detailed balance condition for non-Hamiltonian systems."
   Rep. Math. Phys. 10, 249-258. [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/003448777690046X)
   *Original quantum detailed balance definition.*

2. **Kossakowski, Frigerio, Gorini, Verri (1977)**: "Quantum detailed balance and KMS condition."
   Comm. Math. Phys. 57, 97-110. [Project Euclid](https://projecteuclid.org/journals/communications-in-mathematical-physics/volume-57/issue-2/Quantum-detailed-balance-and-KMS-condition/cmp/1103901281.pdf)
   *KMS-DB connection established.*

3. **Fagnola, Umanita (2007)**: "Generators of Detailed Balance Quantum Markov Semigroups."
   Inf. Dim. Anal. Quantum Prob. 10, 335-363. [arXiv:0707.2147](https://arxiv.org/abs/0707.2147)
   *Modern treatment: GKSL representation of DB generators.*

4. **Buca, Prosen (2012)**: "A note on symmetry reductions of the Lindblad equation."
   New J. Phys. 14, 073007. [arXiv:1203.0943](https://arxiv.org/abs/1203.0943)
   *Weak vs strong Liouvillian symmetry classification.*

5. **Albert, Jiang (2014)**: "Symmetries and conserved quantities in Lindblad master equations."
   Phys. Rev. A 89, 022118. [arXiv:1310.1523](https://arxiv.org/abs/1310.1523)
   *Extended symmetry framework including steady-state structure.*

6. **Medvedyeva, Essler, Prosen (2016)**: "Exact Bethe ansatz spectrum of a tight-binding chain with dephasing noise."
   Phys. Rev. Lett. 117, 137202. [arXiv:1606.09122](https://arxiv.org/abs/1606.09122)
   *η-pairing symmetry; 1D free-fermion ancestor of Π.*

7. **Alhambra, Woods (2017)**: "Dynamical maps, quantum detailed balance, and the Petz recovery map."
   [arXiv:1609.07496](https://arxiv.org/abs/1609.07496)
   *For their QDB dissipative semigroup with no unitary part, the map equals its
   Petz recovery map; an additional commuting unitary part is reversed while
   the same dissipative evolution is retained.*

8. **Roberts, Lingenfelter, Clerk (2021)**: "Hidden Time-Reversal Symmetry, Quantum Detailed Balance and Exact Solutions of Driven-Dissipative Quantum Systems."
   PRX Quantum 2, 020336. [arXiv:2011.02148](https://arxiv.org/abs/2011.02148)
   *Hidden TRS in Lindbladians; closest framework to Π.*

9. **Sá, Ribeiro, Prosen (2023)**: "Symmetry Classification of Many-Body Lindbladians: Tenfold Way and Beyond."
   Phys. Rev. X 13, 031019. [arXiv:2212.00474](https://arxiv.org/abs/2212.00474)
   *38-fold classification of the shifted generator; local-dephasing examples resolve to BDI/CI by parity sector. Π supplies a sectorwise P generator after resolving Π², but our full sector algebra remains uncomputed.*

10. **Kawasaki, Mochizuki, Obuse (2022)**: "Topological phases protected by shifted sublattice symmetry in dissipative quantum systems."
   Phys. Rev. B 106, 035408. [arXiv:2201.09283](https://arxiv.org/abs/2201.09283)
   *Shifted sublattice symmetry: ordinary chiral/sublattice symmetry plus a constant decay-rate shift, eigenvalue pairs ±λ + ia. This is a catalogued comparison for the SHAPE of Π's shifted spectrum. The repository separately establishes an interacting k-body Π and a 6 → 4 → 2 per-site-product-vs-entangled locality ceiling; equivalence and prior-art coverage remain open.*

11. **Haga et al. (2023)**: Incoherenton paper.
    *XY-weight grading = incoherenton number. Bands but no palindrome.*

12. **Chen, Kastoryano, Gilyén (2025)**: "Efficient Quantum Gibbs Samplers with KMS Detailed Balance Condition."
    Comm. Math. Phys. [arXiv:2404.05998](https://arxiv.org/abs/2404.05998)
   *Strict KMS-symmetric Lindbladians for Gibbs sampling; real spectra for that
   symmetric generator.*

13. **Roberts et al. (2025)**: "Hidden Time Reversal in Driven Spin Chains."
    Phys. Rev. Lett. 134, 130404. [PRL](https://doi.org/10.1103/PhysRevLett.134.130404)
    *Extension of hidden TRS to spin chains; new dissipative phase transitions.*

---

## Suggested Next Steps

1. **Formalize the Choi-Jamiolkowski connection.** Map Π from superoperator
   space to an operator in the doubled Hilbert space. Check if it corresponds
   to an antiunitary operator in the Roberts et al. thermofield double framework,
   and do the same for the thermal element M_β of Question 4, whose antilinear
   form K∘Γ_β∘S_XY already carries a bath temperature.

2. **Test Roberts et al. hidden TRS for our system.** Compute whether the
   XXZ + Z-dephasing Liouvillian satisfies their hidden TRS condition, and the
   XXZ + thermal one at one ratio.
   If the first holds, Π and hidden TRS may be dual descriptions of one symmetry;
   if the second, M_β and hidden TRS.

3. **Complete the sectorwise Sá et al. classification.** Work with the centered
   Liouvillian L_c, resolve Π² and every other commuting unitary/strong symmetry,
   then determine the antiunitary generators and their algebra in each irreducible
   block.
   Direct raw-multiset consecutive-gap ratios retain zero gaps and count undefined
   0/0 separately; this degenerate unresolved population has no standard-ensemble
   calibration for random-matrix universality. These statistics do not remove the symmetry class or prove integrability. See
   [Random Matrix Theory](../experiments/RANDOM_MATRIX_THEORY.md).
   The operator result is exact: Π is linear, Π⁴ = I,
   Π² = U_X, and it anticommutes with both L_c and L_c†. After fixing
   U_X parity it phase-normalizes to an involutive P generator. This alone
   does not choose a global or sectorwise SRP class. The gain-loss system
   (Σγ = 0) retains the λ ↔ −λ pairing through its threshold, with two qubits
   per chain generically an exceptional point.
   See [PT-Symmetry Analysis](../experiments/PT_SYMMETRY_ANALYSIS.md).

4. **Contact Roberts/Clerk group and Prosen group.** Both would immediately
   recognize the significance of a shifted anti-similarity. Roberts for
   the hidden TRS connection, Prosen for the η-pairing generalization.
   The classification side of this step has its vehicle: the outbound
   adapter
   [Shifted Order-4 Chiral Symmetry](outbound/SHIFTED_ORDER4_CHIRAL_SYMMETRY.md)
   packages the placement walk of this document, the generalized-P
   resolution, and the statistics data point in the classification
   community's own language.

---

*The Π anti-similarity is not a detailed-balance or KMS condition. It is an
algebraic mirror of the pure-Z dephasing generator; no bath temperature or
thermodynamic equilibrium is implied by Π's mirror; the thermal element M_β of
Question 4 is where a temperature enters. The pure-Z generator may
independently satisfy Alicki/standard QDB with respect to `I/d` after its
Hamiltonian derivation is split off.*
