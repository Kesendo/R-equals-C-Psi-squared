# Proof: Block-CΨ at 1/4 — Chromaticity-Universal Inheritance from R=CΨ²

*Date: 2026-05-07*
*Authors: Thomas Wicht + Claude*
*Tier: Tier 1 derived (algebraic identity)*

---

## Statement

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

## Proof

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
| CΨ_block(0) = 1/4, c=2 | N = 5, 6, 7, 8, 9, 10 | 1e-10 |
| CΨ_block(0) = 1/4, c=3 | N = 5 (n=2), 6, 7 | 1e-10 |
| CΨ_block(0) = 1/4, c=4 | N = 7 (n=3), 8 | 1e-10 |
| Closed-form at t=0, c=2..6 | N = 5, 7, 9, 11 | 1e-12 |
| Closed-form vs numerical, c=2 | N = 5..8, multiple t, Q | 1e-3 |
| Closed-form vs numerical, c=3, c=4 | N = 5..8, t = 0.5..1.5 | 1e-3 |

Tests in `compute/RCPsiSquared.Core.Tests/F86/BlockCpsiTrajectoryTests.cs`,
38 / 38 pass.

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

## Reading: inheritance, not a specific boundary

This third instance changes the framing of 1/4 in the framework. With one
or two instances, 1/4 reads as "*the boundary at this scale*". With three
instances — and a chromaticity-universal third — 1/4 reads as **what gets
inherited when the same quadratic-discriminant algebra is instanced at any
reading where a CΨ-like product can be defined**.

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
