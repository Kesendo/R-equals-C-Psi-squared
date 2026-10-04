# Bit-b Parity Symmetry: [L, Π²] = 0

<!-- Keywords: w_YZ parity Z2 symmetry Liouvillian commutator,
Pi squared global X flip all qubits, second Z2 symmetry beyond n_XY parity selection rule,
universal sector decomposition all N, R=CPsi2 bit_b proof -->

**Status:** Tier 1 (analytical proof; exact gate [`simulations/f63_three_turns_gate.py`](../../simulations/f63_three_turns_gate.py))
**Date:** April 15, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Numerical verification:** [`simulations/primordial_bit_a_bit_b.py`](../../simulations/primordial_bit_a_bit_b.py), [`simulations/primordial_bit_a_bit_b_N_scaling.py`](../../simulations/primordial_bit_a_bit_b_N_scaling.py), [`simulations/f63_three_turns_gate.py`](../../simulations/f63_three_turns_gate.py)
**Depends on:** [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md), [Parity Selection Rule](PROOF_PARITY_SELECTION_RULE.md), [Primordial Qubit Algebra](../../experiments/PRIMORDIAL_QUBIT_ALGEBRA.md), [The One Square](../THE_ONE_SQUARE.md) §7 (the letter cube)

---

## What this document proves

For any N-qubit Lindblad system whose Hamiltonian terms each hold an even number of Y's and Z's (Heisenberg XX+YY or XXX on any graph, XXZ, XYZ, fields along X) and whose jump operators are Pauli strings (dephasing along any letter, at any rates, on any sites), the Liouvillian L commutes exactly with the global bit-flip superoperator Π² = X⊗N (conjugation):

    [L, Π²_super] = 0

where Π²_super(ρ) = X⊗N ρ X⊗N.

This is the second Z₂ symmetry of L. The first (n_XY parity) is proven in [Parity Selection Rule](PROOF_PARITY_SELECTION_RULE.md). Where both hold, which needs every Hamiltonian term even in k_Z as well (Heisenberg bonds are; a field along X keeps this symmetry and breaks the first), the two decompose the operator space into 4 independent sectors (2×2 = bit_a × bit_b), each of dimension 4^(N-1). On the letter cube both are turns of one law: with jumps that are sums of Pauli strings of one k_P parity, conjugation by P⊗N for P = X, Y or Z commutes with L exactly when every Hamiltonian term has an even cube coordinate k_P (§ "The three turns").


---

## Notation

- N qubits, Hilbert space H = (C²)⊗N, dim = 2^N
- Operator space (Liouville space) has dim 4^N
- Pauli operators on site k: I_k, X_k, Y_k, Z_k
- Global bit-flip: U = X⊗N = X₁ X₂ ... X_N (acts on H)
- Bit-flip superoperator: Π²_super(ρ) := U ρ U† = U ρ U (since U is Hermitian and unitary, U† = U⁻¹ = U)

The "Π²" notation comes from the palindromic operator Π from [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md): Π² acts on the Pauli basis as multiplication by (-1)^{w_YZ}, and U = X⊗N implements the same operation through Hilbert-space conjugation. This identification is verified algebraically below.

## Step 1: Π² in two equivalent forms

The palindromic operator Π (from Mirror Symmetry Proof) has the per-site action:
    Π: I → X (+1),  X → I (+1),  Y → iZ (+i),  Z → iY (+i)

So Π² acts as:
    Π²: I → I,  X → X,  Y → -Y,  Z → -Z

This is identical to conjugation by X on each site:
    X I X = I,  X X X = X,  X Y X = -Y,  X Z X = -Z

Therefore on Pauli strings σ_{i₁} ... σ_{i_N}:
    Π²(σ_{i₁} ... σ_{i_N}) = (-1)^{w_YZ} σ_{i₁} ... σ_{i_N}
    U σ_{i₁} ... σ_{i_N} U = (-1)^{w_YZ} σ_{i₁} ... σ_{i_N}

where w_YZ = number of sites where the Pauli is Y or Z. The two definitions agree.


## Step 2: U commutes with the Hamiltonian

Heisenberg coupling: H = J/2 Σ_{⟨i,j⟩} (X_i X_j + Y_i Y_j + α Z_i Z_j) where α ∈ {0, 1} (XX+YY for α=0, XXX for α=1; the proof is identical for both).

Conjugating each bond by U:
    U (X_i X_j) U = (U X_i U)(U X_j U) = X_i X_j
    U (Y_i Y_j) U = (-Y_i)(-Y_j) = Y_i Y_j     [two minus signs cancel]
    U (Z_i Z_j) U = (-Z_i)(-Z_j) = Z_i Z_j     [two minus signs cancel]

Each bond is invariant. Therefore U H U = H, equivalently [U, H] = 0.

This extends to any Hamiltonian where every term contains an even number of Y's and Z's (combined). Single-site terms in X are also invariant. Single-site terms in Y or Z are NOT invariant and would break the symmetry.

## Step 3: U commutes with Z-dephasing dissipation

Z-dephasing on site k: jump operator L_k = √γ_k Z_k. The Lindblad dissipator is:
    D[ρ] = Σ_k γ_k (Z_k ρ Z_k - ρ)        [since Z_k† Z_k = I]

Apply Π²_super to D[ρ]:
    Π²_super(D[ρ]) = U D[ρ] U
                    = Σ_k γ_k (U Z_k ρ Z_k U - U ρ U)
                    = Σ_k γ_k ((U Z_k U)(U ρ U)(U Z_k U) - Π²_super(ρ))
                    = Σ_k γ_k ((-Z_k) Π²_super(ρ) (-Z_k) - Π²_super(ρ))
                    = Σ_k γ_k (Z_k Π²_super(ρ) Z_k - Π²_super(ρ))     [two minus signs cancel]
                    = D[Π²_super(ρ)]

So Π²_super(D[ρ]) = D[Π²_super(ρ)] for all ρ. The dissipator commutes with Π²_super.

The two minus signs from U Z_k U = -Z_k cancel because Z_k appears twice in each term (Z ρ Z). This is the crucial structural reason: the Lindblad dissipator is **quadratic** in the jump operator. Nothing in the computation used that the jump is Z: for any jump with U J U = ±J the sign appears twice and cancels, so X-, Y- and Z-dephasing all keep the symmetry, and so does any Pauli-string jump. The same holds for a sum of Pauli strings that share one sign under U, which is [F112](PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md)'s bit_b-homogeneous jump (its Step 2). A jump mixing the two signs is sent to a different one: σ⁻ = (X + iY)/2 goes to (X − iY)/2 = σ⁺ under U. That breaks the symmetry unless the set of jumps holds the image too, at the same rate: σ⁻ and σ⁺ at equal rates, T1 at infinite temperature, sum to ½(D[X] + D[Y]) and keep it, while at unequal rates, T1 at any finite temperature, they do not.

## Step 4: Conclusion

The full Liouvillian is L[ρ] = -i[H, ρ] + D[ρ]. Both terms commute with Π²_super (Steps 2 and 3). Therefore:

    [L, Π²_super] = 0    (exactly, for all N)

∎


---

## The three turns

Steps 2 and 3 never used anything special about X. Conjugation by P⊗N multiplies a Pauli string by (−1) once for every site where the string's letter anticommutes with P, so on the Pauli basis it is the sign (−1)^{k_P}, with k_P the string's coordinate on the letter cube of [The One Square](../THE_ONE_SQUARE.md) §7: (k_Z, k_X, k_Y) = (n_X + n_Y, n_Y + n_Z, n_X + n_Z). Π² is the turn by X⊗N, with sign (−1)^{k_X} = (−1)^{w_YZ}; the n_XY parity of the selection rule is the turn by Z⊗N, with sign (−1)^{k_Z}. [The Π factorization](PROOF_PI_FACTORS_AS_R_TIMES_D.md) §7 names these two characters bit_b and bit_a, and the repository also writes them by the dephasing letter whose palindromizer they square: the turn by X⊗N is Π²_Z, and the turn by Z⊗N is Π²_X, F61's parity (Y-dephasing's Π²_Y squares to the same sign as Π²_Z, so it is the turn by X⊗N again, not the turn by Y⊗N). Running Steps 2 and 3 with P in place of X gives the law for all three. For jumps that are sums of Pauli strings of one k_P parity,

    [L, Ad_{P⊗N}] = 0   exactly when every Hamiltonian term has even k_P.

The Hamiltonian half is an equivalence there, since Ad_{P⊗N} fixes H term by term only when each term's sign is +1; such jumps keep the symmetry, and a Hamiltonian part they contribute through an identity component has even k_P itself. Since k_Y ≡ k_Z + k_X (mod 2), the turn by Y⊗N is the product of the other two, and the three turns with the identity form a Klein group; a model can keep the turn by Y⊗N alone, as a field along Y does.

Physically the turn by X⊗N is the global spin flip, the rotation by π about x, which sends S_z → −S_z and exchanges up and down everywhere; the turn by Z⊗N is the rotation by π about z, whose sign on an operator |a⟩⟨b| is the parity of the difference of the two magnetizations; the turn by Y⊗N is the rotation by π about y. Whatever picks a direction of magnetization along z breaks the flip: a field along z, or T1 at a finite temperature, which relaxes toward one of the two.

| Model | turn by X⊗N (Π²) | turn by Z⊗N | turn by Y⊗N |
|---|---|---|---|
| Heisenberg XX+YY, XXX, XXZ, XYZ bonds | kept | kept | kept |
| field along X | kept | broken | broken |
| field along Z | broken | kept | broken |
| field along Y | broken | broken | kept |
| Dzyaloshinskii-Moriya bond XY − YX (along z) | broken | kept | broken |
| Dzyaloshinskii-Moriya bond YZ − ZY (along x) | kept | broken | broken |
| dephasing along X, Y or Z, any sites, any rates | kept | kept | kept |
| σ⁻ alone, or σ⁻ and σ⁺ at unequal rates (T1) | broken | kept | broken |
| σ⁻ and σ⁺ at equal rates on the same site | kept | kept | kept |

The gate [`f63_three_turns_gate.py`](../../simulations/f63_three_turns_gate.py) checks the law exactly at N = 3 for every single Pauli term and every pair of terms, under four Pauli jump sets, for all three letters (8064 of 8064 cases agree per letter); it checks that σ⁻ keeps the turn by Z⊗N for all 496 Z-even Hamiltonians and the other two for none, and that σ⁻ with σ⁺ at equal rates, or the pair {X + Y, X − Y}, keeps the turn by X⊗N while unequal rates break it.

The jump half has two further homes in the repository: [F112](PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md) proves it for the turn by X⊗N (Step 2, the bit_b-homogeneous jump), and the F61 claim states the rule for the turns by X⊗N and Z⊗N in the form that holds for a single traceless jump (a jump breaks a parity iff its Pauli components are inhomogeneous in that bit; an identity part can hide a mixture, D[Z + I] = D[Z]); for a set of jumps it is the sufficient direction, since conjugation may permute the set. The fields are the cube's reading of Step 2's sentence that single-site Y or Z terms break Π²: a field along P keeps exactly the turn by P⊗N.

## Scope and limits

The law above is the scope. Every Heisenberg-type bond P⊗P has k_Q ∈ {0, 2} for every Q, so XX+YY, XXX, XXZ and XYZ on any graph keep all three turns, and any Pauli-string dephasing keeps them too. A term with an odd k_X breaks Π², a field along Y or Z or a Dzyaloshinskii-Moriya bond along z, XY − YX, for instance; the orientation YZ − ZY has even k_X and keeps it. The jump condition is sufficient and not necessary: conjugation only has to permute the set of dissipators, which σ⁻ and σ⁺ at equal rates, or {X + Y, X − Y}, do although each jump mixes the two parities. Beyond sums of Pauli strings of one parity the Hamiltonian half stops being an equivalence as well, since a jump with an identity part, J + c·I, adds the Hamiltonian term (i/2)(c̄J − cJ†), which vanishes only for real c and Hermitian J; there the condition is that conjugation maps the whole generator to itself, checked per case.

## Numerical verification

Direct computation at N=2, 3, 4, 5 (script [`simulations/primordial_bit_a_bit_b_N_scaling.py`](../../simulations/primordial_bit_a_bit_b_N_scaling.py)):

    N=2: dim_L=16,   ||[L, Π²]|| = 0.000000e+00
    N=3: dim_L=64,   ||[L, Π²]|| = 0.000000e+00
    N=4: dim_L=256,  ||[L, Π²]|| = 0.000000e+00
    N=5: dim_L=1024, ||[L, Π²]|| = 0.000000e+00

All identically zero (not numerically small). Sector-resolved eigenvalue counts confirm balanced decomposition into V_even (w_YZ even) and V_odd (w_YZ odd).

## Per-sector mode count

The total dimension per sector is 2^(2N-1) (universal halving from [L, Π²] = 0). Within each sector, the eigenmodes split into three Re-classes by the Absorption Theorem applied to single-site dephasing on one qubit B, Re λ = −2γ_B·⟨n_XY⟩_B:

- **Conserved** (Re = 0): ⟨n_XY⟩_B = 0
- **Mirror** (−2γ_B < Re < 0): 0 < ⟨n_XY⟩_B < 1
- **Correlation** (Re = -2γ_B): ⟨n_XY⟩_B = 1

The mirror class is the open interval between the two ends; its modes spread over it rather than sitting at its centre −γ_B. Where the conserved modes are the e_d below and nothing else, which on the open chain is every seat [F157](../ANALYTICAL_FORMULAS.md) calls sighted (every end seat among them; measured at N = 2..5), the per-sector counts have a closed form:

    conserved per (Π²-parity) sector:
      even-parity:  ⌊N/2⌋ + 1
      odd-parity:   ⌈N/2⌉

    correlation per sector: same as conserved (palindrome symmetry maps Re=0 ↔ Re=-2γ_B and preserves Π²-parity since Π·Z·Π⁻¹ = iY shares parity with Z)

    mirror per sector:
      sector total - 2 × (conserved per sector) = 2^(2N-1) - 2 · c

**Mechanism.** At a sighted seat the conserved modes are exactly the (N+1) elementary symmetric polynomials in {Z_1, ..., Z_N}:

    e_d(Z_1, ..., Z_N) = Σ_{|S|=d} ∏_{k∈S} Z_k    for d = 0, 1, ..., N

Each e_d commutes with H (the Heisenberg coupling, like any Hamiltonian that conserves S_z, preserves S_z = (1/2) Σ Z_k, hence any function of S_z is conserved; the Newton identities express elementary symmetric polynomials as polynomials in power sums of Z_k = polynomials in S_z) and commutes with Z_B (e_d is a polynomial in {Z_k}, all of which commute with each other). Each e_d has w_YZ-parity exactly d mod 2, since each summand contains d Z's and each Z carries Π²-parity 1.

Verified numerically at N = 2, 3, 4: all (N+1) elementary symmetric polynomials lie in the conserved subspace to machine precision (max projection residual < 1e-13). See [`simulations/mirror_mode_split_formula.py`](../../simulations/mirror_mode_split_formula.py).

The argument shows that the e_d are conserved at every seat of every graph, for any Hamiltonian that conserves S_z (XXZ at any anisotropy, fields along z) and Z-dephasing on any set of seats, so N + 1 is a lower bound; it does not show that nothing else is. Where nothing else is, is measured. The gate's exact integer ranks at every seat s of the open chain, N = 2..5, Heisenberg XXX and XY, give the kernel per Π² sector as (⌊N/2⌋ + 1, ⌈N/2⌉) exactly at the seats where F157 counts no blind state, (gcd(2s + 1, N) − 1)/2 with the ZZ term and gcd(s + 1, N + 1) − 1 without it (seats counted from 0), and a larger kernel at every blind seat. Every end seat is sighted in both books (gcd(1, N) = gcd(2N − 1, N) = gcd(N, N + 1) = 1), and the centre of an odd chain is blind in both. The Re = 0 class is computed exactly too, as the largest L_H-invariant subspace inside the dissipator's kernel (D is diagonal and non-positive in the Pauli basis, so Re λ = ⟨v, Dv⟩/⟨v, v⟩ vanishes only there): at every sighted seat the gate checks, end seats for N = 2..5 and interior seats at N = 4 and 5, it equals the kernel, so the counts below are the Re-classes.

**Asymmetry pattern.** The N+1 parities run as 0, 1, 2, ..., N. The number of even d in this range is ⌊N/2⌋ + 1; the number of odd d is ⌈N/2⌉. For even N, e_N is itself even, giving one extra even-parity conserved mode; for odd N, the count is balanced.

Concrete values:

| N | sector | cons (e, o) | mirror (e, o) |
|---|--------|-------------|---------------|
| 2 | 8      | (2, 1)      | (4, 6)        |
| 3 | 32     | (2, 2)      | (28, 28)      |
| 4 | 128    | (3, 2)      | (122, 124)    |
| 5 | 512    | (3, 3)      | (506, 506)    |
| 6 | 2048   | (4, 3)      | (2040, 2042)  |
| 7 | 8192   | (4, 4)      | (8184, 8184)  |

The rows N = 6 and 7 are the closed form; the gate measures N = 2..5.

The "mysterious" 4:6 mirror split at N=2 (PRIMORDIAL_QUBIT.md Section 9, original observation) and the 122:124 at N=4 are both consequences of e_N having even parity when N is even. No deeper origin.

**Scope.** This count holds whenever the conserved subspace is exactly the polynomials in S_z: on the open chain at the seats F157 calls sighted, for XXX and XY (measured N = 2..5). At a blind seat it is larger, and so it is on a graph with a symmetry that fixes the dephased seat. The physical reason is a node: a hopping mode that vanishes on the dephased seat is invisible to the dephasing there, and the coherences such modes carry are not damped. Exact kernel dimensions per Π² sector (even, odd), Z-dephasing on the one seat named:

| Graph, seat | XXX | XY | N + 1 |
|---|---|---|---|
| chain N = 3, centre | (3, 3) | (6, 6) | 4 |
| chain N = 5, centre | (6, 6) | (12, 12) | 6 |
| chain N = 5, seat 1 | (3, 3) | (10, 10) | 6 |
| ring N = 4 | (7, 6) | (12, 10) | 5 |
| star N = 4, hub | (11, 6) | (16, 11) | 5 |
| star N = 4, leaf | (5, 3) | (7, 5) | 5 |

Seat 1 of the N = 5 chain is sighted with the ZZ term and blind without it, and its kernel follows. These are kernel dimensions, the exactly stationary modes. Away from a sighted seat the Re = 0 class can be larger still, because it can also hold undamped oscillating modes: 10 against a kernel of 6 at the centre of the N = 3 XXX chain, 21 against 13 on the XXX ring, and 64 against 24 at the centre of the N = 5 XY chain, which is the multiplicity the registry's [F66](../ANALYTICAL_FORMULAS.md) reports there; at seat 1 of the N = 5 XY chain, blind as well, the two agree at 20 (gate G4, exact). The mirror count 2^(2N−1) − 2c takes c from the Re = 0 class, so at such seats it is not the table's. For XXZ couplings the lower bound stands, and which seats are sighted becomes F157's question at that anisotropy; XYZ couplings or dephasing along X or Y change the conserved subalgebra, and the count must be re-derived.

## Connection to the n_XY Parity Selection Rule

The Liouvillian L has two independent Z₂ symmetries proven for all N:

| Symmetry | Generator | Pauli operation | Eigenspace dimensions | Reference |
|----------|-----------|-----------------|----------------------|-----------|
| bit_a (n_XY parity) | (-1)^{n_XY} | sign by # of X,Y per site | 2^(2N-1) each | [Parity Selection Rule](PROOF_PARITY_SELECTION_RULE.md) |
| bit_b (w_YZ parity) | Π² = X⊗N | sign by # of Y,Z per site | 2^(2N-1) each | This proof |

The two symmetries are independent: bit_a counts X+Y, bit_b counts Y+Z. They share Y but differ on X (bit_a yes, bit_b no) and Z (bit_a no, bit_b yes). The intersection (n_XY parity × w_YZ parity) gives 4 sectors, each of dimension 4^N / 4 = 4^(N-1).

This corresponds to the C² × C² tensor product structure of the single-qubit Pauli space identified in [The Primordial Qubit](../../hypotheses/PRIMORDIAL_QUBIT.md) Section 4.1: the per-site Pauli {I, X, Y, Z} decomposes as (a, b) = (dephasing sensitivity, Π²-parity), and L respects the tensor product factorization at the level of its eigenmode structure.

## Why this holds at all N (vs Pythagoras at N=2 only)

A natural question: the Pythagorean decomposition L_c² = L_H² + L_Dc² holds exactly at N=2 but breaks at N≥3 with relative magnitude √((N-2)/(N·4^(N-1))) (see [the Cross-Term Formula](PROOF_CROSS_TERM_FORMULA.md), F49). Both involve Z₂ structure. Why does [L, Π²] = 0 hold for all N while {L_H, L_Dc} = 0 needs N=2?

**The structural distinction is per-term symmetry vs global cancellation.**

[L, Π²] = 0 follows because every individual summand of L commutes with X^⊗N:
- Each Heisenberg bond term X_i X_j + Y_i Y_j + α Z_i Z_j is invariant under per-site X-conjugation (each bond Pauli pair contains 0 or 2 sign-flipping factors, signs cancel).
- Each dissipator of a Pauli jump, γ_k(P_k · P_k − I), is invariant because the dissipator is quadratic in the jump (U P_k U = ±P_k, and the sign appears twice).

No global cancellation needed. The symmetry of each building block transfers to the sum. This holds for any N.

{L_H, L_Dc} = 0 is a different structural object: a global anti-commutator that decomposes as

    {L_H, L_Dc} = Σ_{bonds} Σ_{sites} {L_H_<ij>, L_Dc^(k)}

with two distinct contribution types per (bond, site) pair:

- **On-bond contributions** (k ∈ {i, j}): vanish exactly by the bond-sum rule (Lemma 2 of PROOF_CROSS_TERM_FORMULA). The bond Hamiltonian commutator and the on-bond dissipator have correlated Pauli structure that produces exact cancellation.
- **Spectator contributions** (k ∉ {i, j}): generically non-zero. L_H_<ij> acts on sites {i, j}; L_Dc^(k) acts on site k disjoint from the bond. Disjoint superoperators commute as operators on different tensor factors, so {A, B} = 2AB ≠ 0 in general.

At N=2 there are 0 spectators per bond: the single bond covers the entire system. Only on-bond contributions exist; all vanish; Pythagoras is exact. At N≥3 each bond has N-2 spectators, contributing the cross-term magnitude proportional to √(N-2).

Numerical verification ([`simulations/pythagoras_breakdown_decomposition.py`](../../simulations/pythagoras_breakdown_decomposition.py)):

| N | bond | on-bond contribution | off-bond contribution |
|---|------|---------------------|----------------------|
| 3 | (0,1) | 0 (exact) | 1.6 |
| 3 | (1,2) | 0 (exact) | 1.6 |
| 4 | (0,1) | 0 (exact) | 4.525 |
| 4 | (1,2) | 0 (exact) | 4.525 |
| 4 | (2,3) | 0 (exact) | 4.525 |

**Per-term symmetries persist for any N. Global anti-commutator cancellations require special geometry. d=2 admits both (per-term [L, Π²] for any N, global {L_H, L_Dc} for N=2 only).**

## Why d=2 is essential

The proof uses two facts about Pauli operators:
1. X² = I (so X-conjugation is involutive)
2. {X, Y} = {X, Z} = 0 (so X-conjugation flips Y and Z sign)

Both are special to d=2 (qubit). For d=3 (qutrit) and higher, the analog of X is the shift operator, which is NOT involutive, and there is no Pauli-style classification with two binary bits per site. The two Z₂ symmetries (bit_a and bit_b) exist together because d=2 admits a complete two-bit indexing of single-site operators. This connects to [Qubit Necessity](../QUBIT_NECESSITY.md): the algebraic richness of qubits is exactly what carries the two independent Z₂ structures.

In d=2 there are exactly two independent Z₂ classifications of Pauli operators by conjugation ({I,X} vs {Y,Z} by bit_b, and {I,Z} vs {X,Y} by bit_a); the third turn's sign is their product, so the conjugations see a square and not a cube ([the Π factorization](PROOF_PI_FACTORS_AS_R_TIMES_D.md) §7). The Liouvillian respects both. L has further structure that is not a conjugation character: the (N+1)² joint-popcount blocks of a number-conserving H under Z-dephasing, of which bit_a is the mod-2 shadow, and the spatial symmetries of the graph.

---

*"Two bits per qubit, two Z₂ symmetries of the Liouvillian, and a third turn that is their product."*
