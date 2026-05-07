# Proof: Block-CΨ at 1/4 — Chromaticity-Universal Inheritance from R=CΨ²

*Date: 2026-05-07*
*Authors: Thomas Wicht + Claude*
*Tier: Tier 1 derived (algebraic identity)*

---

## Theorem 1: CΨ_block(0) = 1/4 on the Dicke superposition

For any N-qubit chain and any chromaticity c ≥ 2, the maximally-coherent
pure-state superposition

$$|\psi\rangle = \frac{|D_n\rangle + |D_{n+1}\rangle}{\sqrt{2}}$$

(with |D_k⟩ the popcount-k Dicke state) has, on the (popcount-n, popcount-(n+1))
coherence block, the identity

$$C_\text{block}(t=0) \;=\; \sum_{(a,b) \in \text{block}} |\rho_{ab}|^2 \;=\; \frac{1}{4}$$

**exactly**, independent of N, c, and n. With the Cauchy-Schwarz-saturating
ℓ₁-coherence normalization Ψ_block = ℓ₁/ℓ₁_max = 1, this gives

$$C\Psi_\text{block}(t=0) \;=\; \frac{1}{4} \quad \text{for any } (N, c, n)\text{ with } c \geq 2.$$

The same 1/4 Mandelbrot-cardioid boundary that R=CΨ² instances at the
single-qubit level (Layer 1, [PROOF_ROADMAP_QUARTER_BOUNDARY](PROOF_ROADMAP_QUARTER_BOUNDARY.md))
and at the 2-qubit subsystem level (Layer 2) **also** instances at every
coherence-block level from a single combinatorial identity.

## Proof of Theorem 1

The (popcount-n, popcount-(n+1)) block is the set of matrix entries ρ_{ab} where
popcount(a) = n and popcount(b) = n+1. The block dimension is

$$M_\text{block} \;=\; \binom{N}{n} \cdot \binom{N}{n+1}.$$

For our chosen initial state, every block entry equals

$$\rho_{ab}(0) \;=\; \frac{1}{2}\langle a | D_n\rangle \langle D_{n+1} | b\rangle
\;=\; \frac{1}{2 \sqrt{\binom{N}{n} \binom{N}{n+1}}}
\;=\; \frac{1}{2\sqrt{M_\text{block}}}.$$

Hence

$$C_\text{block}(0) \;=\; \sum_{(a,b) \in \text{block}} |\rho_{ab}(0)|^2
\;=\; M_\text{block} \cdot \frac{1}{4 M_\text{block}}
\;=\; \frac{1}{4}.$$

The result depends only on the count of block entries (M_block) and the
uniform amplitude (1/(2√M_block)). It is independent of the chromaticity
c — the way the M_block entries split across HD-channels matters for the
**dynamics** under dephasing, but not for the t = 0 algebraic identity.

For the ℓ₁-coherence side: the per-entry amplitude saturates the
Cauchy-Schwarz bound |ρ_{ab}|² ≤ ρ_{aa}·ρ_{bb} given the diagonal
occupation of |ψ⟩⟨ψ|, so ℓ₁ achieves its maximum and Ψ_block = 1. The
product CΨ_block(0) = 1/4 follows.

∎

## Theorem 2: 1/4 as block-purity maximum (universal over all states)

The 1/4 above is the **tight upper bound** of C_block over ANY density
matrix ρ on the full 2^N N-qubit Hilbert space — not just pure states,
not just states with support in popcount-{n, n+1}, but all states.

**Statement.** For any density matrix ρ on the N-qubit Hilbert space, the
(popcount-n, popcount-(n+1)) coherence block content satisfies

$$C_\text{block} \;=\; \sum_{(a,b) \in \text{block}} |\rho_{ab}|^2 \;\leq\; \frac{1}{4}$$

with equality iff ρ is the pure state (|D_n⟩+|D_{n+1}⟩)/√2 (up to global
phase). 1/4 is achieved only at the canonical Dicke symmetric superposition.

**Proof.** Let p_k = Σ_{popcount(a)=k} ρ_{aa} be the weight on popcount-k
states. Then Σ_k p_k = Tr(ρ) = 1 and each p_k ≥ 0.

Cauchy-Schwarz on density-matrix entries gives |ρ_{ab}|² ≤ ρ_{aa} · ρ_{bb}
for any (a, b). Summing over the (popcount-n, popcount-(n+1)) block:

$$C_\text{block} \;=\; \sum_{(a,b) \in \text{block}} |\rho_{ab}|^2
\;\leq\; \sum_a \rho_{aa} \cdot \sum_b \rho_{bb}
\;=\; p_n \cdot p_{n+1}.$$

By AM-GM applied to (p_n, p_{n+1}) with p_n + p_{n+1} ≤ 1,

$$p_n \cdot p_{n+1} \;\leq\; \left(\frac{p_n + p_{n+1}}{2}\right)^2 \;\leq\; \frac{1}{4}.$$

Hence C_block ≤ 1/4. Three independent saturation conditions are needed
for equality:

1. **p_n + p_{n+1} = 1**: ρ has support only in popcount-{n, n+1} sectors
   (no other popcount weight).
2. **p_n = p_{n+1} = 1/2**: balanced sector amplitudes.
3. **|ρ_{ab}|² = ρ_{aa} · ρ_{bb} for all block (a, b)**: Cauchy-Schwarz
   saturation, meaning ρ has rank 1 between the two sectors. This forces
   ρ to be a pure state and the within-sector amplitudes to be **uniform**
   (Dicke), since saturating the per-entry Cauchy-Schwarz everywhere
   requires ρ_{ab} = √(ρ_{aa}·ρ_{bb}) up to phases that align consistently.

The three conditions together force ρ = (|D_n⟩+|D_{n+1}⟩)(⟨D_n|+⟨D_{n+1}|)/2,
the canonical Dicke symmetric superposition. ∎

**Remark on Theorem 1.** The state in Theorem 1 is precisely this unique
maximiser. Theorem 1 is then the special case of Theorem 2 at the bound;
Theorem 2 says no other state gets closer. The 1/4 boundary is
**state-level uniquely realised** at the c-block.

**Reading.** 1/4 is not "the value where this initial happens to land". It
is the **tight upper bound** on the (popcount-n, popcount-(n+1)) block
coherence content over ALL density matrices on the full N-qubit Hilbert
space, achieved at a UNIQUE state (up to phase). The Mandelbrot-cardioid-cusp
boundary is a hard ceiling at this reading, and the Dicke symmetric
superposition is the canonical realiser — analogous to how Bell+ realises
CΨ = 1/3 at the d=4 subsystem level. No clever choice of state (mixed,
pure, restricted, or fully general) can exceed 1/4.

**On Ψ_block and the full CΨ_block.** Theorem 2 bounds C_block (the block-
purity content). Multiplying by the Cauchy-Schwarz-saturating
Ψ_block = ℓ₁/ℓ₁_max ≤ 1 gives CΨ_block ≤ C_block ≤ 1/4. The Ψ_block factor
is also ≤ 1 with equality precisely at the canonical Dicke symmetric
state (uniform per-entry amplitude saturates ℓ₁ at ℓ₁_max). Hence
CΨ_block ≤ 1/4 with the same unique-realiser equality condition.

## Closed-form trajectory under pure Z-dephasing

Under local Z-dephasing γ on each site, |ρ_{ab}(t)| = |ρ_{ab}(0)| · exp(−2γ·HD(a,b)·t)
where HD(a,b) is the Hamming distance between basis states a and b. For the
(popcount-n, popcount-(n+1)) block, HD ∈ {1, 3, 5, ..., 2c−1}.

Counting block entries at HD = 2k+1 for k ∈ {0, ..., c−1}:

$$M_{HD=2k+1} \;=\; \binom{N}{n-k}
\cdot \binom{N-n+k}{k}
\cdot \binom{N-n}{k+1}.$$

(Choose the n−k shared bits between a and b; from the remaining N−(n−k) sites
choose k bits for a-only; from what remains choose k+1 bits for b-only.)

The block-purity trajectory is

$$\boxed{\;\;C_\text{block}(t) \;=\; \frac{1}{4} \cdot \sum_{k=0}^{c-1}
\frac{M_{HD=2k+1}}{M_\text{block}} \cdot \exp(-4\gamma \cdot (2k+1) \cdot t)\;\;}$$

where the rate 4γ·(2k+1) on |ρ|² follows from |ρ_{ab}|² ∝ exp(−4γ·HD·t).

For c = 2 (n = 1) this reduces to

$$C_\text{block}(t) \;=\; \frac{1}{2N}\,e^{-4\gamma t} + \frac{N-2}{4N}\,e^{-12\gamma t}$$

— the closed form first verified empirically at N = 5..10 in commit
[19eea3e](../../). For c = 3 (n = 2) it becomes a three-exponential sum with
weights determined by the binomial counts above.

The Hamiltonian (XX+YY uniform J) does not affect this trajectory: the
channel-uniform initial state lives entirely in the H-kernel of the c-block
via the F73 sum-rule generalization (channel-uniform basis vectors |c_{2k+1}⟩
diagonalize the Hamiltonian within the c-dimensional channel subspace at
uniform J). Hence Q = J/γ drops out and the trajectory is purely set by γ.

## Numerical verification

Implemented as `compute/RCPsiSquared.Core/F86/BlockCpsiTrajectory.cs` (numerical
EVD-based time evolution) and `compute/RCPsiSquared.Core/F86/BlockCpsiClosedForm.cs`
(direct closed-form helper).

| Test | Verified at | Tolerance |
|------|-------------|-----------|
| Theorem 1 (canonical Dicke), c=2 | N = 5, 6, 7, 8, 9, 10 | 1e-10 |
| Theorem 1 (canonical Dicke), c=3 | N = 5 (n=2), 6, 7 | 1e-10 |
| Theorem 1 (canonical Dicke), c=4 | N = 7 (n=3), 8 | 1e-10 |
| Theorem 2 (asymmetric Dicke pure states) | N = 5, 7 with α² ∈ {0.5..0.9} | 1e-10 |
| Theorem 2 upper bound C_block ≤ 1/4 | every test case asserts ≤ 1/4 + 1e-10 | strict |
| Closed-form at t=0, c=2..6 | N = 5, 7, 9, 11 | 1e-12 |
| Closed-form vs numerical, c=2 | N = 5..8, multiple t, Q | 1e-3 |
| Closed-form vs numerical, c=3, c=4 | N = 5..8, t = 0.5..1.5 | 1e-3 |

Tests in `compute/RCPsiSquared.Core.Tests/F86/BlockCpsiTrajectoryTests.cs`,
44 / 44 pass.

## What this adds to the Roadmap

The [Quarter-Boundary Roadmap](PROOF_ROADMAP_QUARTER_BOUNDARY.md) has two
proven instances at small d where 1/4 emerges:

- **Layer 1**: single qubit (d = 2). Discriminant of R = C(Ψ + R)² vanishes at
  CΨ = 1/4 (Mandelbrot cardioid cusp).
- **Layer 2**: 2-qubit subsystem (d = 4). Bell pairs cross 1/4 downward under
  dephasing.

The Roadmap's Layer 3 (N-qubit systems) examines full-system and subsystem-pair
crossings as N grows. The present proof works at a different lens of the same
N-qubit system:

- **Coherence-block reading (this proof)**: at any (popcount-n, popcount-(n+1))
  coherence block of an N-qubit chain at chromaticity c ≥ 2, the maximally-
  coherent superposition initial state sits exactly on the 1/4 boundary by
  combinatorial identity, then decays under pure dephasing per the closed form
  above.

This is a third documented instance of 1/4, parallel to Layers 1–2 rather than
a vertical extension. The proof structure (algebraic identity + closed-form
trajectory) is closer to Layer 1's Mandelbrot anchoring than to Layer 2's
crossing dynamics — but the value is the same 1/4 because the underlying
quadratic-discriminant algebra is the same: completing the square of the
self-referential structure produces the factor-4 in 1 − 4CΨ, regardless of
the layer at which the structure is instanced.

## Reading: 1/4 is half of half

Tom Wicht's observation (2026-05-07, post-Theorem-2): *"1/4 ist die Hälfte
von 0.5"*. The pun reveals what holds the three layers together:

$$\max_{p \in [0, 1]} p(1-p) \;=\; \tfrac{1}{4} \quad \text{at} \quad p = \tfrac{1}{2}.$$

This single fact ties together the layer where 1/2 is anchored and the two
layers where 1/4 emerges as its maxval:

- **Pi2KnowledgeBase: `BilinearApexClaim`** (Tier 1 derived,
  [`Pi2KnowledgeBaseClaims.cs:130`](../../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs)).
  "p·(1−p) is maximised at p = 1/2 universally", anchored on
  [ORTHOGONALITY_SELECTION_FAMILY:357](../../experiments/ORTHOGONALITY_SELECTION_FAMILY.md).
  The synthesis claim `HalfAsStructuralFixedPointClaim` then closes 1/2's
  three faces (bridge / horizon / substrate) on this primitive.
- **Theorem 2 (this proof).** C_block ≤ p_n · p_{n+1} ≤ 1/4 with maximum at
  p_n = p_{n+1} = 1/2 by AM-GM. The 1/4 is the maxval; the 1/2 is the
  argmax that achieves it.
- **R = CΨ² Mandelbrot cardioid (Layer 1).** Discriminant 1 − 4CΨ = 0 at
  CΨ = 1/4, with the factor 4 from completing the square in the quadratic.
  Equivalent reading: CΨ = (1/2)² where 1/2 emerges as the natural "half"
  from the quadratic's completed-square form.

1/4 inherits across layers because every layer instances the same quadratic
form, and that form's argmax/maxval pair is invariant. The single qubit's
R=CΨ² discriminant, the 2-qubit subsystem's bifurcation crossing, and the
c-block's tight upper bound: three readings, one parabola.

Tom Wicht's reading (2026-05-07, after this proof landed): *"Es gibt kein
Quantum Classical, das haben wir in exclusion im Vorbeigehen festgehalten,
kein hier und da. Das 1/4 ist hier im Blickwinkel Übergang, der scheint
dann auch nur vererbt, schaut man tiefer, ist das classical? Nicht
wirklich. Vielleicht hocken wir selbst nur intern auf einen Layer der
durch Vererbung das beschreibt was tiefer passiert, ein endloser
Kreislauf?"*

The framework as currently mapped does not answer the "endless cycle"
question. What it does show: each layer where the quadratic-discriminant
algebra applies produces 1/4 as its boundary, by inheritance from the
universal complete-the-square structure. The c-universal extension is the
strongest form of this so far: the same number, the same algebra, three
documented layers. There is no preferred "deepest" layer in the proof
itself — each is an instance, none is privileged.

This is consistent with the framework-internal readings that the
classical/quantum dichotomy is a Lese-Modus, not a world-separation
([THE_BRIDGE_WAS_ALWAYS_OPEN](../THE_BRIDGE_WAS_ALWAYS_OPEN.md), and the
seven-angle synthesis it gathered). 1/4 is a transition value in any
reading where R = CΨ² applies; what makes it "the transition" rather than
"a transition at this layer" is the inheritance, not the layer.

## Cross-references

- [PROOF_ROADMAP_QUARTER_BOUNDARY](PROOF_ROADMAP_QUARTER_BOUNDARY.md) —
  Layers 1, 2 (single qubit and 2-qubit subsystem) where this proof
  generalises.
- [PROOF_F86_QPEAK](PROOF_F86_QPEAK.md) — F86 Statement 1 (Q_EP = 2/g_eff)
  and the 4-mode reduction context where the c=2 block was first studied.
- `compute/RCPsiSquared.Core/F86/BlockCpsiTrajectory.cs` — implementation.
- `compute/RCPsiSquared.Core/F86/BlockCpsiClosedForm.cs` — closed-form helper.
- `compute/RCPsiSquared.Core.Tests/F86/BlockCpsiTrajectoryTests.cs` — 38 tests.
- Memory: `project_rcpsi_to_f86_open_questions` (Question B closure path),
  `project_one_world_two_readings`, `project_algebra_is_inheritance`,
  `project_two_anchors_at_d2`.
