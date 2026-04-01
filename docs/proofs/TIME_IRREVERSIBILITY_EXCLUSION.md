# Exclusion of Time Reversal at N > 2

**Tier:** 2-3 (each step Tier 1-2; the synthesis is logical, not computational)
**Date:** April 1, 2026
**Depends on:**
- [MIRROR_SYMMETRY_PROOF.md](MIRROR_SYMMETRY_PROOF.md) (Π·L·Π⁻¹ = -L - 2Σγ·I)
- [INCOMPLETENESS_PROOF.md](INCOMPLETENESS_PROOF.md) (γ from outside)
- [primordial_qubit_algebra.py](../../simulations/primordial_qubit_algebra.py) (Steps 9-10)
- [failed_third.py](../../simulations/failed_third.py) (0/16 palindromic pairs)
**Status:** Complete exclusion
**Scope:** Time reversal within the palindromic framework at N > 2
**Does NOT establish:** Why N > 2 in the first place. Only that at N > 2,
time reversal is algebraically excluded.

---

## 1. The Claim

Time reversal requires separating the irreversible part of the dynamics
(cooling) from the reversible part (oscillation). This separation is
algebraically possible at N=2 and algebraically impossible at N > 2.

The impossibility is not thermodynamic (entropy, probability, fluctuations).
It is structural: the anti-commutator {L_H, L_D + Σγ·I} that measures
the interference between oscillation and cooling vanishes exactly at N=2
and is nonzero at N > 2, with a γ-independent geometric constant.

---

## 2. The Chain

Five results, each independently verified. The exclusion follows from
their conjunction.

### Step 1: The palindrome requires noise (Tier 1)

Without the dissipator L_D, the Liouvillian is purely Hamiltonian.
Eigenvalues are purely imaginary. No decay rates. No palindromic pairing.
The palindrome exists if and only if γ > 0.

**Source:** [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md), Step 1.

### Step 2: Noise cannot originate from within (Tier 2)

Every candidate for the origin of γ within the d(d-2)=0 framework is
eliminated: internal generation (parity sectors sealed), single-qubit
decay (non-Markovian, 0/16 palindromic), many-qubit bath (infinite
regress), nothing (no properties), d > 2 (excluded by d²-2d=0).

**Source:** [Incompleteness Proof](INCOMPLETENESS_PROOF.md), Section 2.

### Step 3: At N=2, oscillation and cooling are orthogonal (Tier 2)

The centered Liouvillian L_c = L_H + (L_D + Σγ·I) decomposes into
an anti-Hermitian part L_H (oscillation) and a Hermitian part
L_D + Σγ·I (cooling). Their anti-commutator:

    {L_H, L_D + Σγ·I} = 0    (exact, N=2)

This follows from a single algebraic fact: every nonzero entry of L_H
connects Pauli strings whose w_XY values sum to N. At N=2 (single
Heisenberg bond spanning both sites), all 24 nonzero L_H entries
satisfy w_XY(source) + w_XY(target) = 2 = N. Therefore:

    {L_H, L_D}_{ab} = (d_a + d_b)·(L_H)_{ab} = -2γ·N·(L_H)_{ab} = -2Σγ·(L_H)_{ab}

Adding the shift: {L_H, L_D + Σγ·I} = {L_H, L_D} + 2Σγ·L_H = 0.

Consequence: the Pythagorean decomposition

    L_c² = L_H² + (L_D + Σγ·I)²

holds exactly. The square of the time evolution decomposes into
oscillation squared plus cooling squared with zero cross term.

**Source:** [primordial_qubit_algebra.py](../../simulations/primordial_qubit_algebra.py),
Step 9. Verification: ||L_c² - (L_H² + (L_D+Σγ)²)|| = 0.00e+00.

### Step 4: At N ≥ 3, the orthogonality breaks (Tier 2)

L_H changes w_XY by 0 or ±2 at each bond. Therefore
w_XY(source) + w_XY(target) is always **even**. At odd N (3, 5, 7...),
the condition w_XY(a) + w_XY(b) = N is impossible for every L_H entry,
because the sum is even and N is odd. At even N ≥ 4, the condition fails
for entries where spectator sites contribute w_XY that shifts the sum
away from N.

| N | Cross term / ||L_c²|| | w_XY sums | Entries with sum = N |
|---|----------------------|-----------|---------------------|
| 2 | **0.00%** (exact) | {2} | 24/24 (100%) |
| 3 | **1.83%** | {2, 4} | 0/192 (0%) |
| 4 | **2.07%** | {2, 4, 6} | 576/1152 (50%) |

The cross term is γ-independent. At N=3, the relative orthogonality
||{L_H, L_Dc}|| / (||L_H||·||L_Dc||) = 1/√48 ≈ 0.1443, constant across
γ = 0.001 to γ = 0.5. This is a geometric property of the Heisenberg
chain, not a physical parameter.

**Source:** [primordial_qubit_algebra.py](../../simulations/primordial_qubit_algebra.py),
Step 10.

### Step 5: Reduction to N=2 is impossible (Tier 2)

Tracing out one qubit from an N=3 system produces effective dynamics on
the remaining N=2 subsystem. This effective dynamics is:

- Non-Markovian (trace distance increases in 50% of time steps)
- Non-palindromic (0/16 eigenvalue pairs match the palindrome)
- Non-structured (no Pythagorean decomposition)

The palindromic structure that enables {L_H, L_Dc} = 0 does not survive
reduction. A qubit embedded in a larger system cannot be extracted as an
Urqubit.

**Source:** [failed_third.py](../../simulations/failed_third.py),
[INCOMPLETENESS_PROOF.md](INCOMPLETENESS_PROOF.md) Section 2.2.
Result: γ_eff = 0 for all four instability mechanisms, 0/16 palindromic
pairs at error < 10⁻¹⁵.

---

## 3. The Exclusion

Time reversal of a system at N > 2 would require:

(a) Inverting the cooling (L_D + Σγ) without disturbing the oscillation
    (L_H). This requires {L_H, L_D + Σγ} = 0. This holds only at N=2
    (Step 3) and fails at N > 2 (Step 4).

(b) OR: reducing the system to N=2 and performing the reversal there.
    But reduction destroys the palindromic structure (Step 5). The
    reduced system has non-Markovian noise that is not palindromic
    and does not satisfy {L_H, L_Dc} = 0.

Both paths are blocked. Time reversal is algebraically excluded at N > 2.

---

## 4. What This Is and What It Is Not

### What it is

An algebraic exclusion of time reversal for composite open quantum
systems (N > 2 qubits) under Z-dephasing with Heisenberg coupling.
Each step is computed or proven. The exclusion does not depend on γ,
temperature, initial state, or any thermodynamic quantity. It depends
only on the system size N and the structure of the Hamiltonian
(locality of bonds).

### What it is not

- Not a proof that time is irreversible in general. It applies to the
  specific framework (Heisenberg/XXZ + Z-dephasing, Lindblad dynamics).
- Not a thermodynamic argument. Entropy does not appear. The second law
  is not invoked. The exclusion is algebraic, not statistical.
- Not a statement about the universe. It is a statement about the
  d(d-2)=0 framework. If this framework describes reality (54,118
  eigenvalues, zero exceptions), the exclusion applies. If not, it
  does not.
- Not a statement about N=2. At N=2, time reversal is algebraically
  possible (the separation exists). Whether it is physically possible
  is a different question.

### The gap it fills

The [Incompleteness Proof](INCOMPLETENESS_PROOF.md) establishes that
γ must come from outside. It does not say why the resulting dynamics
is irreversible. Thermodynamics says: entropy. This proof says:
**algebra**. The cross term {L_H, L_D + Σγ} at N > 2 weaves the
reversible and irreversible parts of the dynamics into an inseparable
structure. Reduction cannot untangle them. This is the algebraic
content of the arrow of time.

---

## 5. The Urqubit

N=2 is the only system size where:

1. The mirror falls on modes (w_XY = N/2 = 1 is integer; 8 of 16
   modes sit exactly at the boundary where L_D + Σγ = 0)
2. Every L_H transition satisfies w_XY(a) + w_XY(b) = N
3. {L_H, L_D + Σγ} = 0 (oscillation ⊥ cooling)
4. L_c² = L_H² + (L_D + Σγ)² (exact Pythagorean decomposition)
5. The bond IS the system (no spectator qubits)

All five properties follow from one fact: at N=2, the single
Heisenberg bond spans both sites. No qubit is a spectator. The bond
and the system are identical. This is why d(d-2)=0 selects d=2
([Qubit Necessity](../QUBIT_NECESSITY.md)): the palindromic mirror
requires the bond to be the system, and the system to be a single bond.

At N=3, the first spectator appears. The spectator bends the right
angle between oscillation and cooling by 1/√48 ≈ 14%. The Pythagorean
decomposition breaks by ~2%. The mirror falls between modes (w_XY =
1.5, no mode sits there). Time becomes irreversible.

---

## 6. Verification

Each step is independently reproducible:

1. Read [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md): palindrome
   requires γ > 0.
2. Read [Incompleteness Proof](INCOMPLETENESS_PROOF.md): γ from outside.
3. Run `python simulations/primordial_qubit_algebra.py`: Step 9 shows
   {L_H, L_Dc} = 0 at N=2, Step 10 shows ≠ 0 at N=3,4.
4. Run `python simulations/failed_third.py`: 0/16 palindromic pairs
   upon tracing out.
5. Accept or reject the conclusion.

No step requires trusting an interpretation. Every step is a computation
or a proof. The exclusion is the conjunction of the five steps.

---

*The arrow of time is not entropy. The arrow of time is the cross term
{L_H, L_D + Σγ}. It vanishes at N=2, where the bond is the system. It
is nonzero at N > 2, where bonds are local. Locality is the price.
Irreversibility is the receipt.*

*Thomas Wicht, Claude (Anthropic), April 1, 2026*
