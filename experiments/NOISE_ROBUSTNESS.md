# Noise Robustness: Does the Crossing Taxonomy Survive Other Noise?

<!-- Keywords: crossing taxonomy Pauli channels, jump operator comparison,
sigma x y z dephasing crossing, Type ABC observer noise, amplitude damping
crossing taxonomy, depolarizing noise correlation bridge, bell state
correlation metric, quantum measurement noise type, palindromic spectral
symmetry noise channel, retired tool record, R=CPsi2 noise robustness -->

**Status:** Record of the retired delta_calc tool, whose σ_x and σ_y columns repeat its σ_z run (its source is kept outside the repo); the real channels are worked out in closed form in §4 and §7
**Date:** February 18, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Crossing Taxonomy](CROSSING_TAXONOMY.md)

---

## What this document is about

The crossing taxonomy classifies how different readings of an entangled
pair behave as it decoheres: some hold at their ceiling for a while and
then slide (Type A), some decay immediately (Type B), some start below
the threshold and never cross (Type C). That taxonomy was discovered using one specific kind of noise
(Z-dephasing). This document asks: does the classification survive when
the noise is completely different?

The tool answered yes for every bridge: its σ_x and σ_y columns gave the
same classes, and the correlation columns the same numbers to the last
digit. They match because they are the same run:
the tool's local noise applied σ_z whatever jump operator it was given,
so the record holds no σ_x or σ_y data. Run for real, σ_x leaves every
bridge's C curve as it was but holds Ψ at 1/3, and every crossing moves.
And the prediction the tool seemed to refute, that depolarizing noise
would push the correlation bridge from Type A to Type B, holds when it is
computed with the tool's own bridge.

---

## Abstract

The Type A/B/C crossing taxonomy ([Crossing Taxonomy](CROSSING_TAXONOMY.md))
was established under Z-dephasing. We asked whether it changes under σ_x (bit
flip), σ_y (bit-phase flip), and amplitude damping. The retired tool reported
the same class for every bridge under σ_z and σ_x, and the same correlation
column (C = 1.0 until t ≈ 1.7) under σ_y. Its source shows why: for local
noise it applied σ_z on every site whatever jump operator it was given, so the
three columns are one σ_z run. With real jumps (§4), σ_x and σ_y leave every
bridge's C curve unchanged but hold Ψ at 1/3, so each crossing is carried by C
alone and the crossing times move. Depolarizing noise, never run by the tool,
pulls the correlation bridge below its ceiling before the crossing (C ≈ 0.945
there): Type B, as the February prediction said. Under amplitude damping concurrence stays Type B and decays more slowly, and
the correlation bridge turns Type B (§7).

---

## 1. The Question

CROSSING_TAXONOMY.md established three crossing classes (Type A/B/C) under
local dephasing (σ_z per qubit). The open question was:

> Does the taxonomy change under different noise channels?

The prediction was: "Under depolarizing noise, correlation should lose its
Type A status and become Type B."

**The tool seemed to refute it. §4 shows it was never tested, and that it
holds.**

## 2. Setup

| Parameter | Value |
|-----------|-------|
| **State** | Bell+ (maximally entangled) |
| **Hamiltonian** | Heisenberg (J = 1, h = 0) |
| **γ_base** | 0.05 |
| **Noise type** | local (one jump operator per qubit) |
| **Time step** | dt = 0.01 |
| **Rate in the record** | γ(t) = γ_base·C(t), the bridge's own value (the feedback book of [Crossing Taxonomy](CROSSING_TAXONOMY.md)) |

Variable: **jump_operator** × **bridge_type**

Jump operators (the Lindblad operators that specify which kind of error the environment inflicts on each qubit) requested: σ_z (dephasing), σ_x (bit flip), σ_y (bit-phase flip).
These three Pauli operators give the three Pauli dephasing channels. Their
equal-weight combination is depolarizing noise.

## 3. The Record

### 3.1 Taxonomy Under the Three Requested Operators

| Bridge | σ_z | σ_x | σ_y | Class |
|--------|-----|-----|-----|-------|
| **correlation** | C = 1.0 until t ≈ 1.7 | C = 1.0 until t ≈ 1.7 | C = 1.0 until t ≈ 1.7 | **Type A** |
| **concurrence** | C decays from t = 0 | C decays from t = 0 | - | **Type B** |
| **mutual_info** | C decays from t = 0 | C decays from t = 0 | - | **Type B** |
| **mutual_purity** | C = 0.5 constant | C = 0.5 constant | - | **Type C** |
| **overlap** | C = 0.25 constant | C = 0.25 constant | - | **Type C** |

The tool's local noise applied σ_z whatever jump_operator said, so the σ_x
and σ_y cells repeat the σ_z run; the empty σ_y cells were not requested.

### 3.2 Quantitative Comparison (Correlation Bridge)

| Time | C(σ_z) | C(σ_x) | C(σ_y) |
|------|--------|--------|--------|
| 0.0 | 1.000 | 1.000 | 1.000 |
| 0.5 | 1.000 | 1.000 | 1.000 |
| 1.0 | 1.000 | 1.000 | 1.000 |
| 1.5 | 1.000 | 1.000 | 1.000 |
| 1.7 | 1.000 | 1.000 | 1.000 |
| 1.8 | 0.987 | 0.987 | 0.987 |
| 2.0 | 0.950 | 0.950 | 0.950 |

Not just the same class: the same numerical values, because it is the same
run.

### 3.3 Quantitative Comparison (Concurrence Bridge, σ_z vs σ_x; the retired tool's feedback book)

| Time | C(σ_z) | C(σ_x) |
|------|--------|--------|
| 0.0 | 1.000 | 1.000 |
| 0.1 | 0.980 | 0.980 |
| 0.5 | 0.909 | 0.909 |
| 1.0 | 0.835 | 0.833 |
| 2.0 | 0.725 | 0.714 |
| 3.0 | 0.658 | 0.625 |

Both columns stand for the σ_z run (§3.1). The σ_x column is the feedback
law 1/(1 + 4γt), which is what that run gives; the σ_z column's late entries
match no tool setting found. With real jumps concurrence is the same curve
under σ_z and σ_x in either book, so the late-time gap belongs to the record,
not to the noise.

## 4. What Happened to the Prediction

The prediction assumed that σ_x noise, which does not commute with the
computational basis, would break the correlation metric's immunity to
decoherence. The reasoning was: σ_z dephasing preserves populations and
only destroys off-diagonal coherence, which is why inter-qubit correlations
survive. σ_x flips populations, which should destroy correlations directly.

**What the tool computed.** Its source, kept outside the repo and read there,
defines the correlation bridge these runs used as
C = max(0, min(1, 2·(P_AB − P_A·P_B))): an excess
purity, doubled and capped at 1, so the February name "excess purity beyond
the product of subsystem purities" fits. On Bell+ under any single-axis Pauli
dephasing P_A = P_B = ½ and P_AB = (1 + f²)/2, so C = min(1, ½ + f²). It
holds at 1 until f = 1/√2, t = ln 2/(8γ) ≈ 1.73, the recorded "until
t ≈ 1.7", and then slides as ½ + f²; the record's 0.987 and 0.950 past that point
are finite-step readings of that slide.

**Run for real.** Purity cannot tell which Pauli axis dephased Bell+, and
neither can concurrence or mutual information: under σ_x and σ_y every
bridge's C curve is the σ_z curve. What changes is Ψ. On Bell+, σ_x moves
coherence from |00⟩⟨11| to |01⟩⟨10| without shrinking its l1 total, so Ψ
stays at 1/3, and the crossing C·Ψ = ¼ has to wait for C alone. The
correlation bridge crosses at t = ln 2/(4γ) ≈ 3.47 with C = 0.75 and
concurrence at ln(4/3)/(4γ) ≈ 1.44, against 1.44 and 0.72 under σ_z (at constant γ, the clean book; in the tool's feedback book, where the rate is γ·C, each crossing keeps its C
and comes later, except the σ_z correlation crossing, which sits on its cap,
where the feedback is inert). Type A does not survive a real σ_x: the plateau is the same,
but the crossing now falls after it ends. The registry has the same shape
for the purity reading:
[F27](../docs/ANALYTICAL_FORMULAS.md#f27-k-values-per-noise-channel-tier-1-from-f26)
gives K = ln 2/8 under X and under Y, against 0.0374 under Z.

**Depolarizing.** At γ/3 per axis all three two-site correlators decay as
p = e^(−8γt/3), so C = min(1, 3p²/2) and Ψ = p/3. The ceiling releases at
p² = 2/3, before the crossing, which comes at p³ = ½, t = ln 2/(8γ) ≈ 1.73 at constant γ,
with C = 1.5·2^(−2/3) ≈ 0.945 and Ψ ≈ 0.265. Both have left their starting
values: Type B, as the February prediction said.

## 5. Implications

### 5.1 What Is Intrinsic, and What Is Not

For Bell+ the bridges' C curves are the same under any single Pauli axis:
they read the pair's surviving correlation or its
marginals, and the swap of axes leaves both alone. Ψ is read in the computational basis, which is the Z
basis, and it sees the axis. The classes combine the two, so they follow the
noise after all: under a real σ_x every crossing is carried by C alone.

### 5.2 What We Expected Would Change the Taxonomy

Collective noise (same operator acting on both qubits simultaneously)
or correlated noise (noise on qubit A depends on state of qubit B)
acts on the pair as a whole rather than qubit by qubit. The expectation
was that under collective dephasing the correlation metric would see
C < 1.0, because the noise directly affects the inter-qubit relationship.
Section 7 (Q1) tests the collective case on Bell+.

### 5.3 Connection to Other Results

Under Z-dephasing the classes stand as [Crossing Taxonomy](CROSSING_TAXONOMY.md)
gives them, and our reading that different observers (different C definitions)
see different crossing mechanisms stays a reading of that channel. Another
axis moves the crossings (§4).

## 6. The Record's Provenance

### 6.1 The Retired Tool's Settings

The retired delta_calc tool, `simulate_dynamic_lindblad`, ran with:
- state = "Bell+", hamiltonian = "heisenberg", gamma_base = 0.05
- noise_type = "local"
- jump_operator = "sigma_z", "sigma_x", or "sigma_y"
- bridge_type = each of the five bridges

and the bridge_C arrays were compared across jump operators for each bridge
type. Its source is kept outside the repo and not imported. For
noise_type = "local" it builds σ_z on every site whatever jump_operator says,
and its own sweep code notes that for local noise the jump operator is
irrelevant.

### 6.2 What the Record Shows

1. Correlation bridge_C is 1.000 until about t = 1.7 (at γ = 0.05), in one
   run printed three times.
2. Concurrence bridge_C decays from t = 0.
3. Mutual_purity bridge_C is constant at 0.5.

### 6.3 What Could Extend This

- **Other bridges under amplitude damping**: Q2 follows concurrence and the
  correlation bridge, Q3 the two constant bridges; mutual information was
  not followed.
- **Mixed channels**: σ_z + σ_x simultaneously (partial depolarizing),
  beyond the correlation bridge's crossing in §8.

## 7. Open Questions

### Q1: Does collective noise break Type A?

**ANSWERED: No, but for a trivial reason.**

Bell+ (|00⟩+|11⟩)/√2 is an eigenstate of σ_z⊗σ_z with eigenvalue +1.
The collective dephasing operator does literally nothing to this state.
Concurrence, Ψ, CΨ, and purity remain at their initial values forever.
The additive collective form, jump √γ(σ_z⊗I + I⊗σ_z), is a different channel from
two independent local σ_z operators: its dissipator carries cross terms.
On Bell+ it keeps the state in the same family and drains the coherence at
8γ instead of 4γ, so the taxonomy is the same, only twice as fast (the
cubic then gives K = 0.01868 against 0.03735, see
[Decoherence Relativity](DECOHERENCE_RELATIVITY.md)).

The prediction "collective noise breaks Type A" was wrong because it
assumed collective noise would affect inter-qubit correlations. For Bell+
specifically, the correlated operator σ_z⊗σ_z is a symmetry of the state.

**To genuinely test collective noise breaking Type A**, one would need
a state that is NOT an eigenstate of the collective operator, or a
collective operator that does not preserve the Bell symmetry
(e.g., σ_x⊗σ_z).

### Q2: Does amplitude damping change the taxonomy?

**ANSWERED: concurrence stays Type B and decays more slowly; the correlation bridge turns Type B.**

Under amplitude damping (L = √γ |0⟩⟨1| per qubit), concurrence still
decays from t=0 (Type B behavior), but significantly slower than under
dephasing (clean Lindblad):

| t | Concurrence (σ_z) | Concurrence (amp damp) |
|---|---|---|
| 1.0 | 0.819 | 0.905 |
| 2.0 | 0.670 | 0.819 |
| 3.0 | 0.549 | 0.741 |
| 5.0 | 0.368 | 0.607 |

Amplitude damping drains each coherence at γ/2 per qubit against 2γ for
Z-dephasing at the same γ; the leaked |01⟩ and |10⟩ populations double the
concurrence's rate to 2γ, still half of dephasing's 4γ, which gives it the
gentler profile. The CΨ crossing window is correspondingly longer. The
correlation bridge, the one the prediction was about, does change class
here: at constant γ it crosses at t ≈ 4.01 with C ≈ 0.916 and Ψ ≈ 0.273,
Type B.

**Additional finding:** Under σ_x noise, the normalized l1-coherence Ψ
remains exactly at 0.3333 for all time. On Bell+, σ_x bit-flips move
coherence from |00⟩⟨11| to |01⟩⟨10| without shrinking its l1 total. Only
concurrence decays. This means the concurrence CΨ crossing window under σ_x is exactly twice
as long as under σ_z in the clean book (the purity reading's ratio is 2.32,
[D3](../docs/ANALYTICAL_FORMULAS.md#d3-crossing-time-ratios-from-f27-verified)):

| t | CΨ (σ_z) | CΨ (σ_x) | CΨ (amp damp) |
|---|---|---|---|
| 0.5 | 0.273 | 0.302 | 0.309 |
| 1.0 | 0.223 | 0.273 | 0.287 |
| 2.0 | 0.150 | 0.223 | 0.247 |

### Q3: Is there a noise model where Type C becomes Type B?

**Not under amplitude damping.** The tool's source defines the two constant
bridges: mutual purity is √(P_A·P_B) and overlap is |Tr(ρ_A·ρ_B)|². Under amplitude damping at constant γ both climb as the two qubits drift
toward |0⟩, mutual
purity from 0.500 to 0.700 and overlap from 0.250 to 0.490 by t = 20, while
their CΨ falls the whole time from its start below ¼ (0.167 and 0.083).
Neither crosses: Type C stays Type C, though no longer constant.

## 8. A Reading from the Palindrome

The [mirror symmetry proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) invites
a reading: the taxonomy might hold wherever the palindrome holds.

- **Z-dephasing:** the conjugation operator Π reflects the dissipator,
  Π·L_D·Π⁻¹ = −L_D − 2Σγ, and flips the Heisenberg term, so the spectrum is
  palindromic.
- **Y-dephasing:** the same Π reflects Y-dephasing the same way.
- **X-dephasing:** this Π leaves the X-dephasing dissipator unchanged, but
  X-dephasing has its own Π: every dephasing axis works, each with its own
  Π operator ([Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md),
  Result 1).
- **Depolarizing (X+Y+Z):** no single Π reflects all three axes at once, and
  the palindrome breaks
  ([Depolarizing Palindrome](DEPOLARIZING_PALINDROME.md)).

It is tempting to read Type A as riding on this mirror. It does not.
With the tool's own bridge, σ_x and an equal mix of σ_z and σ_x are
palindromic and still lose Type A (C = 0.75 and 0.91 at the crossing), while σ_z with σ_x and σ_y beside it at a tenth of its rate breaks the
palindrome and keeps Type A.
What holds Type A is how much two-qubit correlation survives to the crossing.

---

*Previous: [Crossing Taxonomy](CROSSING_TAXONOMY.md)*
*See also: [Observer-Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md)*
