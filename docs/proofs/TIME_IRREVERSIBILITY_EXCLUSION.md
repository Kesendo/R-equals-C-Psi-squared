# Loss of Oscillation-Cooling Orthogonality at N > 2 (and the Arrow-of-Time Reading)

**Tier:** the orthogonality result (Steps 1–5) is Tier 1–2 (computed/proven); the
"time-reversal exclusion / arrow of time" READING is **Tier 3 interpretation** (see the §3 caveat).
**Date:** April 1, 2026 (reframed 2026-06-22, see §3 caveat)
**Depends on:**
- [MIRROR_SYMMETRY_PROOF.md](MIRROR_SYMMETRY_PROOF.md) (Π·L·Π⁻¹ = -L - 2Σγ·I)
- [INCOMPLETENESS_PROOF.md](INCOMPLETENESS_PROOF.md) (γ from outside)
- [primordial_qubit_algebra.py](../../simulations/primordial_qubit_algebra.py) (Steps 9-10)
- [failed_third.py](../../simulations/failed_third.py) (0/16 palindromic pairs)
- [PROOF_CROSS_TERM_FORMULA.md](PROOF_CROSS_TERM_FORMULA.md) (R(N) = √((N−2)/(N·4^(N−1))))
**Status:** the orthogonality / Pythagorean-decomposition result is complete and exact.
**Scope:** the **Frobenius-Pythagorean orthogonality** of the oscillatory (L_H) and cooling
(L_D+Σγ·I) parts of the Liouvillian: their anti-commutator vanishes exactly at N=2 and is lost
(by a γ- and topology-independent constant) for N>2.
**Does NOT establish:** that time reversal is *literally* excluded. The DYNAMICAL separability
obstruction (whether the flow e^{(L_H+L_Dc)t} factors into oscillation × cooling) is the
**commutator** [L_H, L_Dc] (BCH: e^{A+B}=e^A e^B iff [A,B]=0), which is **nonzero at *all* N,
including N=2** (‖[L_H,L_Dc]‖≈22.6 at N=2). The vanishing of the *anti*-commutator at N=2 is a
Frobenius-orthogonality fact, not a separability/reversibility criterion, so the "arrow of time"
reading below is an interpretive identification, not a time-reversal theorem.

---

## What this document is about

At N=2 qubits, the oscillatory and dissipative parts of the dynamics are
perfectly orthogonal **in the Frobenius inner product**: the anti-commutator
{L_H, L_D+Σγ·I} vanishes exactly, and L_c² splits as L_H² + (L_D+Σγ)² with zero
cross term. At N > 2, this orthogonality breaks, the cross term is nonzero, with a
fixed geometric magnitude R(N) independent of γ and topology. That is the rigorous,
self-contained content of this document.

We READ this N-dependent loss as the algebraic root of an arrow of time (the cross term
"weaves" oscillation and cooling together once spectators appear). That reading is a
**Tier-3 interpretation**, not a separability theorem: *literal* dynamical separability
(undoing cooling without disturbing oscillation) is governed by the **commutator**
[L_H, L_Dc], which is nonzero at *every* N, including N=2, so the flows never truly
factor, even where the anti-commutator vanishes. The honest claim is about Frobenius
orthogonality and its loss, not about whether time can be reversed.

---

## 1. The Claim

Time reversal requires separating the irreversible part of the dynamics
(cooling) from the reversible part (oscillation). This separation is
algebraically possible at N=2 and algebraically impossible at N > 2.

The impossibility is not thermodynamic (entropy, probability, fluctuations).
It is structural: the anti-commutator (the symmetrized product
{A,B} = AB + BA, which vanishes when two operators are "orthogonal")
{L_H, L_D + Σγ·I} that measures the interference between oscillation
and cooling vanishes exactly at N=2
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

| N | Cross term / \|\|L_c²\|\| | w_XY sums | Entries with sum = N |
|---|----------------------|-----------|---------------------|
| 2 | **0.00%** (exact) | {2} | 24/24 (100%) |
| 3 | **1.83%** | {2, 4} | 0/192 (0%) |
| 4 | **2.07%** | {2, 4, 6} | 576/1152 (50%) |

The cross term is γ-independent. At N=3, the relative orthogonality
||{L_H, L_Dc}|| / (||L_H||·||L_Dc||) = 1/√48, following the general
formula R(N) = √((N-2)/(N·4^(N-1))) proven in
[PROOF_CROSS_TERM_FORMULA](PROOF_CROSS_TERM_FORMULA.md). The constant
is topology-independent (chain = ring = star = complete) and depends
only on N.

**Source:** [PROOF_CROSS_TERM_FORMULA](PROOF_CROSS_TERM_FORMULA.md),
[primordial_qubit_algebra.py](../../simulations/primordial_qubit_algebra.py) Step 10.

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

## 3. The Reading, and its caveat

The arrow-of-time reading goes: time reversal of a system at N > 2 would require:

(a) Inverting the cooling (L_D + Σγ) without disturbing the oscillation
    (L_H), which the document associates with {L_H, L_D + Σγ} = 0 (holds only
    at N=2, Step 3; fails at N > 2, Step 4).

(b) OR: reducing the system to N=2 and performing the reversal there,
    but reduction destroys the palindromic structure (Step 5).

So the Frobenius orthogonality (a) and the clean reduced structure (b) both hold
**only at N=2**, and one reads this as an arrow of time appearing at N>2.

> **Caveat (the load-bearing one, added 2026-06-22).** Step (a)'s premise, that
> "inverting cooling without disturbing oscillation requires the *anti*-commutator
> {L_H, L_Dc} to vanish", is **not correct as a separability criterion**, and the
> "exclusion" therefore does not follow as a theorem. Whether the flow e^{(L_H+L_Dc)t}
> factors into a pure-oscillation piece and a pure-cooling piece (i.e. whether one can
> undo one without the other) is governed by the **commutator** [L_H, L_Dc] (Baker-
> Campbell-Hausdorff: e^{A+B}=e^A e^B iff [A,B]=0), **not** the anti-commutator. And
> **‖[L_H, L_Dc]‖ ≈ 22.6 ≠ 0 already at N=2**, the very point this document calls
> "separable / reversible". So oscillation and cooling do *not* dynamically separate even
> at N=2; the anti-commutator vanishing there is a Frobenius-orthogonality fact with no
> bearing on reversibility. (Independently: the full Liouvillian spectrum is exactly
> palindromic at N=2, 3, 4 alike, and e^{L_c t} is invertible at every N, §7, so nothing
> distinguishes N=2 from N>2 regarding actual time-reversal.) **What survives** is the real,
> γ- and topology-independent result of Steps 3–4: the oscillation-cooling Frobenius-
> Pythagorean orthogonality holds only at N=2 and is lost for N>2 (constant R(N)). The
> "time reversal excluded at N>2 / the arrow of time is the cross term" framing is a Tier-3
> **interpretation** of that geometric fact, not a derived exclusion.

---

## 4. What This Is and What It Is Not

### What it is

An exact algebraic result on the **Frobenius-Pythagorean orthogonality** of
oscillation and cooling for composite open quantum systems under Z-dephasing with
Heisenberg coupling: the anti-commutator {L_H, L_D+Σγ} vanishes only at N=2 and is lost
for N>2, by a γ- and topology-independent constant R(N) that depends only on N and bond
locality. (Plus a Tier-3 *reading* of that loss as an algebraic arrow of time, see the
§3 caveat for why the reading is an interpretation, not a time-reversal exclusion.)

### What it is not

- Not a proof that time is irreversible in general. It applies to the
  specific framework (Heisenberg/XXZ + Z-dephasing, Lindblad dynamics).
- Not a thermodynamic argument. Entropy does not appear. The second law
  is not invoked. The exclusion is algebraic, not statistical.
- Not a statement about the universe. It is a statement about the
  d(d-2)=0 framework. If this framework describes reality (87,376
  eigenvalues, zero exceptions), the exclusion applies. If not, it
  does not.
- Not a statement that N=2 is reversible. At N=2 the *Frobenius* orthogonality
  {L_H, L_Dc}=0 holds, but the *dynamical* separation still does not exist there
  (the commutator [L_H, L_Dc]≠0), so even N=2 is not "reversible" in the literal sense.
  N=2 is special only for the anti-commutator orthogonality, not for time reversal.

### The gap it fills

The [Incompleteness Proof](INCOMPLETENESS_PROOF.md) establishes that
γ must come from outside. Thermodynamics reads the resulting irreversibility as
entropy. This document offers a different *reading*: the cross term {L_H, L_D + Σγ}
at N > 2 weaves the oscillatory and dissipative parts into a Frobenius-inseparable
structure once spectators appear. **We call that the algebraic content of the arrow of
time, as a Tier-3 interpretation.** (It is not a derived irreversibility: the dynamics
is dissipative, hence irreversible in the ordinary CPTP sense, at *every* N including
N=2; the cross term is a geometric N-marker of *where* oscillation and cooling stop being
Frobenius-orthogonal, not a switch that turns reversibility off at N>2.)

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
angle between oscillation and cooling by R(3) = 1/√48 ≈ 14% (from the
general formula R(N) = √((N-2)/(N·4^(N-1))),
[proof](PROOF_CROSS_TERM_FORMULA.md)). The Pythagorean
decomposition breaks by ~2%. The mirror falls between modes (w_XY =
1.5, no mode sits there). We read this as the onset of an arrow of time
(Tier-3 interpretation, §3 caveat), the *Frobenius* orthogonality is lost;
literal reversibility was never present (the commutator [L_H,L_Dc]≠0 at all N).

---

## 6. Verification

Each step is independently reproducible:

1. Read [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md): palindrome
   requires γ > 0.
2. Read [Incompleteness Proof](INCOMPLETENESS_PROOF.md): γ from outside.
3. Run `python` [`simulations/primordial_qubit_algebra.py`](../../simulations/primordial_qubit_algebra.py): Step 9 shows
   {L_H, L_Dc} = 0 at N=2, Step 10 shows ≠ 0 at N=3,4.
4. Run `python` [`simulations/failed_third.py`](../../simulations/failed_third.py): 0/16 palindromic pairs
   upon tracing out.
5. Accept or reject the conclusion.

No step requires trusting an interpretation. Every step is a computation
or a proof. The exclusion is the conjunction of the five steps.

---

---

## 7. Computational Reversal: the Information Window

Physical time reversal is excluded at N > 2 (Section 3). But
computational reversal (calculating the past from the present) is
mathematically defined:

    ρ(0) = e^{-Lt} ρ(t)

The backward evolution amplifies every mode by e^{d_k t}, where
d_k = -Re(λ_k) is the decay rate. The fastest mode (d_max = 2(N-1)γ)
determines the precision requirement:

    t_max = p · ln(10) / d_max

where p is the number of available decimal digits of precision.

### The information window (computed, N=2-6, γ=0.05, J=1)

| N | d_max | 64-bit (t_max) | 128-bit | t_steady (info gone) |
|---|-------|---------------|---------|---------------------|
| 2 | 0.20 | 184 J⁻¹ | 391 J⁻¹ | 50 J⁻¹ |
| 3 | 0.30 | 123 J⁻¹ | 261 J⁻¹ | 50 J⁻¹ |
| 4 | 0.40 | 92 J⁻¹ | 196 J⁻¹ | 50 J⁻¹ |
| 5 | 0.50 | 74 J⁻¹ | 157 J⁻¹ | 50 J⁻¹ |
| 6 | 0.60 | 61 J⁻¹ | 131 J⁻¹ | 50 J⁻¹ |

t_max: how far back you can compute at the given precision.
t_steady = 5/d_min: when the system reaches steady state and
information is irretrievably gone.

### Physical scales

| System | d_max | 64-bit window | 128-bit window |
|--------|-------|--------------|----------------|
| N=10 qubits, γ=10⁹ s⁻¹ | 2×10¹⁰ s⁻¹ | 2 ns | 4 ns |
| Macroscopic (N~10²³) | 2×10³² s⁻¹ | 10⁻³¹ s | 10⁻³¹ s |

For a macroscopic system, even with 10¹⁰⁰ digits of precision,
the computable window is ~10⁶⁸ seconds for the fastest mode.
The slow modes survive longer: d_min = 2γ gives a window of
~10⁻⁸ s (the T₂ coherence time) per digit of precision.

### The palindrome as forgetting schedule

The palindromic pairing d ↔ 2Σγ - d determines which information
is lost first and which survives longest:

- Fast modes (d near 2(N-1)γ): information lost first. These are the
  high-w_XY modes, the most "quantum" coherences.
- Slow modes (d near 2γ): information lost last. These are the
  low-w_XY modes, near the classical sector.
- The palindromic partner of a fast-dying mode is a slow-surviving
  mode. The two carry complementary information.

The ratio of survival times: t_slow/t_fast = (N-1). At N=4:
the slowest mode's information lives 3× longer than the fastest.

### What this means

Physical reversal: impossible at N > 2 (cross term, algebraic).

Computational reversal: possible in principle, with cost:
- Precision: p digits give t_max = p·ln(10)/d_max
- The palindrome determines which information survives how long
- After steady state: information is genuinely gone, no computation
  recovers it

The information window has been measured experimentally for decades.
It is called the **coherence time T₂**. The framework identifies it
as the maximum computable depth of the past: the time beyond which
no amount of computation can reconstruct what happened.

---

*(Tier-3 reading.) Where thermodynamics reads the arrow of time as entropy, this
framework offers a geometric image: the oscillation-cooling cross term {L_H, L_D + Σγ}
vanishes at N=2, where the bond is the system, and is nonzero at N > 2, where bonds are
local. Locality is the price; the lost Frobenius orthogonality is the receipt. (Read as
an identification, not a derivation, see §3: the dynamics is dissipative, hence
irreversible, at every N; the cross term marks where oscillation and cooling cease to be
Frobenius-orthogonal, not where time stops being reversible.)*

*Thomas Wicht, Claude (Anthropic), April 1, 2026*
