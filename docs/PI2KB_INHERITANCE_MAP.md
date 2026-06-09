# Pi2 Knowledge Base — Inheritance Map

**Question Tom 2026-05-16:** *"Ich denke das durch Vererbung diese Anker vielleicht auf irgendeine Art zusammenhängen, irgendwo müssen sich die Anker vererben."*

**Answer:** yes — the structural anchors of the R=CΨ² framework are organised as a tree in `compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs`. This document is the visual map of that tree plus the remaining prose-only edges that are not yet typed as Claim relationships.

---

## The tree (typed Claim graph)

```
PolynomialFoundationClaim                  Pi2KnowledgeBaseClaims.cs:441
│   d²−2d=0 ↔ R=CΨ² ; minimum-memory equation
│   The trunk that generates both branches.
│
├── number-branch (1/d = 1/2 lineage)
│   │
│   QubitDimensionalAnchorClaim             :18
│   │   1/d = 1/2 at d=2 (the unsigned magnitude basis)
│   │
│   ├── BilinearApexClaim                   :135
│   │   │  1/2 = argmax of p·(1−p)
│   │   │
│   │   ├── QuarterAsBilinearMaxvalClaim    :176
│   │   │      1/4 = (1/2)² = maxval of p·(1−p)
│   │   │
│   │   └── ArgmaxMaxvalPairClaim           :228
│   │          closure: (1/2 argmax, 1/4 maxval) inseparable invariant
│   │
│   ├── HalfAsStructuralFixedPointClaim     :325
│   │      three faces of 1/2 close (bridge, horizon, substrate)
│   │
│   └── Pi2DyadicLadderClaim                Pi2DyadicLadderClaim.cs:?
│       │   a_n = 2^(1−n) ladder rungs
│       │
│       ├── PolynomialDiscriminantAnchorClaim   PolynomialDiscriminantAnchorClaim.cs:64
│       │      a_{−1} = 4 = polynomial discriminant of d²−2d=0
│       │      mirror partner: a_3 = 1/4 (the QuarterAsBilinearMaxval) — closure 4·(1/4)=1
│       │
│       ├── Pi2OperatorSpaceMirrorClaim     Pi2OperatorSpaceMirrorClaim.cs
│       │
│       ├── AbsorptionTheoremClaim          AbsorptionTheoremClaim.cs
│       │      α = 2γ₀·⟨n_XY⟩  (quantum decay rate)
│       │      └── UniversalCarrierClaim    UniversalCarrierClaim.cs
│       │
│       └── 53 × F*Pi2Inheritance leaves    compute/RCPsiSquared.Core/Symmetry/F*Pi2Inheritance.cs
│              F1, F25, F33, F38, F49b, F60, F62, F65, F73, F75, F76, F80,
│              F81, F82, F83, F84, F85, F86 family, F87, F88a, F88b, F89, F90,
│              F91, F92, F93, …
│
├── angle-branch (90° lineage)
│   │
│   NinetyDegreeMirrorMemoryClaim           Pi2KnowledgeBaseClaims.cs:381
│   │   90° rotation back to the mirror (the i in F80's 2i)
│   │
│   ├── Pi2InvolutionClaim                  :73   ←── NEW TYPED EDGE
│   │   │   Π² commutes with L (squared from F1)
│   │   │   ──── ctor injects F1PalindromeIdentity (2026-05-16) ────
│   │   │
│   │   └── KleinFourCellClaim              :108
│   │          F88a two-axis (Π²_Z, Π²_X) → 4 Klein cells
│   │
│   └── Pi2I4MemoryLoopClaim                Pi2I4MemoryLoopClaim.cs
│          i⁴ = 1 (Z₄ closure of 90° rotations)
│
├── polarity branch
│   PolarityLayerOriginClaim                :516
│       +0 / 0 / −0 polarity layer at d=0 → ±0.5 pair at d=2
│
└── perspectival branch  ←── NEW 2026-05-16
    TwoReadingsClaim                        TwoReadingsClaim.cs
        any object on d²=4^N admits two coordinate readings of one
        underlying object. Names seven layers of this pattern:
        number/angle of d, argmax/maxval of p·(1−p), M/Π·M·Π⁻¹ (F81),
        bra/ket of ρ, inside/outside Q = J/γ₀, classical/quantum Lese-Modus,
        inter-sectoral wave ("we are the standing wave between").
```

Both branches **reconverge** at `Pi2DyadicLadderClaim` (rungs hold 1/2 lineage values) and at the F-leaves (each F-claim cites both Dyadic + I4 Memory Loop).

---

## Tier-1 roots

Three roots in the typed graph (no inbound edges from other Claims):

1. **`PolynomialFoundationClaim`** (`Pi2KnowledgeBaseClaims.cs:441`) — semantic root of the number-anchor lineage. Its docstring states: *"This is the trunk that generates both framework anchors: `HalfAsStructuralFixedPointClaim` (1/2 number-anchor) and `NinetyDegreeMirrorMemoryClaim` (90° angle-anchor)."* Only markdown citations.

2. **`F1PalindromeIdentity`** (`compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs:24`) — root of the angle/Π-branch via its squaring step `Π² · L · Π⁻² = L`. Now explicitly injected into `Pi2InvolutionClaim` as a typed ctor parent (2026-05-16); before that the edge was prose only.

3. **`ChiralKClaim`** (`compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs`) — chiral / sublattice symmetry root. K = diag((−1)^ℓ) anti-commutes with any NN tight-binding Hamiltonian; gives spectrum inversion `E_{N+1−k} = −E_k` at the Hamiltonian level. Sibling to PolynomialFoundation (dimensional) and F1 (Liouvillian); none derives from another. Added 2026-05-16 to close the prose-only edge.

---

## Remaining prose-only edges (asserted but not yet typed)

Each line is an edge that EXISTS in narrative form (docstring, proof, reflection) but is NOT yet a typed Claim-to-Claim ctor injection.

| Edge | Asserted at | Why not typed (yet) |
|---|---|---|
| ~~**F71 spatial mirror → F1 / Π cluster**~~ (map originally misnamed parent) | `MIRROR_SYMMETRY_PROOF.md`, `PROOF_C1_MIRROR_SYMMETRY.md` | **PARTIALLY CLOSED 2026-05-16:** F71's actual prose parents are `HalfIntegerMirrorClaim` (N-blocked from ctor injection — F71 is N-universal, HalfIntegerMirror takes N at construction) and `Pi2DyadicLadderClaim` (N-agnostic). The latter is now a typed ctor parent of `F71MirrorSymmetryPi2Inheritance`. The HalfIntegerMirror edge remains prose-only because of the N-mismatch. |
| ~~**K chiral / sublattice (BDI class)**~~ | ~~`ChiralK.cs`, memory `project_chiral_partnership`~~ | **CLOSED 2026-05-16:** `ChiralKClaim` typed as a sibling root (see Tier-1 roots list above) |
| ~~**F86, F57, F64 as siblings of 1/4**~~ (map was imprecise) | ~~`docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md`, `docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md`~~ | **CLOSED 2026-05-16 (Wave 5):** F57DwellTimeQuarterPi2Inheritance and F64CavityModeExposurePi2Inheritance already took `QuarterAsBilinearMaxvalClaim`. The F86 leg closed via `IbmBlockCpsiHardwareTable` (lowest blast-radius F86 entry, single Lazy factory in F86KnowledgeBase.cs:248). Two heavier F86 → 1/4 leaves (BlockCpsiTrajectory, F86TPeakPi2Inheritance mirror-partner) remain deferred for blast-radius reasons but not for prose-only reasons. |
| ~~**"Two readings" anchor**~~ | ~~`reflections/ON_BOTH_SIDES_OF_THE_MIRROR.md`, `reflections/ON_THE_HALF.md`~~ | **CLOSED 2026-05-16:** `TwoReadingsClaim` typed as Tier1Derived child of `PolynomialFoundationClaim`; enumerates seven layers of the pattern (number/angle, argmax/maxval, M/Π·M·Π⁻¹, bra/ket, inside/outside, classical/quantum, inter-sectoral) |
| ~~**F1 → F81 / F80 chain**~~ (F71 not in this chain after closer reading) | ~~`docs/proofs/PROOF_F81_PI_CONJUGATION_OF_M.md`~~ | **PARTIALLY CLOSED 2026-05-16:** F80 + F81 now take `F1PalindromeIdentity` as a typed ctor parent (their Π and M derive from F1's master statement). F82 + F84 inherit F1 transitively through F81. Remaining F-leaves that use sine modes / specific dynamics (F75, F76, F33, F65) do not directly use Π and so do not get F1 as a typed parent. |
| ~~**F60 → PolarityLayerOrigin + Quarter**~~ | ~~F60 docstring, anchor string — both already named these as "registration discards"~~ | **CLOSED 2026-05-16 (Wave 1):** all three previously-discarded parents (PolarityLayerOrigin, QuarterAsBilinearMaxval, ArgmaxMaxvalPair) are now typed ctor parents. F60 is the first F-formula whose primary anchor sits directly on the 0.5-shift axis. |
| ~~**F62 → Quarter + F61BitAParity**~~ | ~~F62 docstring: "Combined with the Parity Selection Rule (F61), this proves that single-excitation states never cross CΨ = 1/4."~~ | **CLOSED 2026-05-16 (Wave 2):** both previously-discarded parents (Quarter, F61) typed. F62 + F61 together close the SE-regime "structurally outside framework" reading. |
| ~~**F25 → Quarter; F57 → F25 mother-claim**~~ | ~~F25 docstring: "F25 is the mother claim of F57" — F57's prefactor 1.080088 = 2/1.851701 = 2/|dCΨ/dt|_{t_cross} from F25~~ | **CLOSED 2026-05-16 (Wave 3):** F25 → Quarter typed; F57 → F25 typed (mother-claim edge mirrors F75 → F77 pattern); F57 → ArgmaxMaxvalPair also typed. |
| ~~**F77 → Half; F38 → Half**~~ | ~~F77 docstring: "1/2 = a_2 = Bloch baseline = HalfAsStructuralFixedPointClaim"; F38 docstring: "HalfHalfBalance = 1/2 = a_2 = HalfAsStructuralFixedPointClaim"~~ | **CLOSED 2026-05-16 (Wave 4):** both take `HalfAsStructuralFixedPointClaim` as typed parent. F38 ctor change cascaded cleanly through F61/F62/F63 chains via co-author null-checks. F77 → F75 (Taylor-limit mother claim) still deferred as registration discard — F75-chain construction (F66 → F65 → F75) is heavy. |

## Remaining open edges (after the 2026-05-16 sweep)

Only two structural items remain:

1. **F71 → HalfIntegerMirror** — structurally blocked by N-mismatch. `F71MirrorSymmetryPi2Inheritance` is N-universal; `HalfIntegerMirrorClaim` takes `int N` at construction. The edge cannot be a ctor parent without either parameterising F71 by N or introducing a per-N adapter. Acceptable as prose-only.

2. **F77 → F75 (Taylor-limit mother)** — held as registration discard. F75-chain construction (`F66 → F65 → F75`) is multi-step; a small builder helper would tame the cost. Promotion deferred for later.

Three additional F86 → 1/4 edges (BlockCpsiTrajectory, F86TPeakPi2Inheritance mirror-partner via `_ladder.Term(3)`, BlockCoherenceContent helper) live in Sammelbecken territory and were intentionally NOT touched in Wave 5; each carries blast-radius (3 public factories, multi-test surface) that the IbmBlockCpsiHardwareTable closure does not.

## Next direction: anchor-side blind-spot sweep

The closures above all start from a claim and ask "what should be its parent?" (top-down). The remaining methodological complement: start from each anchor and ask "where is my constant used without me being a typed parent?" (bottom-up). Each algebraic handshake needs a visible anchor, otherwise no consumer of the graph can find it. Per-anchor reverse search is the natural Wave 6.

---

## How to read the tree

Edges in the typed graph mean: the child Claim's constructor accepts the parent Claim and stores/uses it (e.g. `PolynomialDiscriminantAnchorClaim` takes three parents in its ctor). When you instantiate a child, you must pass its parents — Pi2KnowledgeBase's constructor (`Pi2KnowledgeBase.cs:73`) does this wiring centrally.

For F-leaves: each `F*Pi2Inheritance.cs` class has a ctor that asks for `Pi2DyadicLadderClaim` (and often `Pi2I4MemoryLoopClaim`, `HalfAsStructuralFixedPointClaim`, …). The wiring is in the `compute/RCPsiSquared.Runtime/PolarityArchitecture/` registration files.

---

## The single numerical structure underneath

All anchors collapse to **one quadratic parabola** `p · (1 − p)`:

- argmax at `p = 1/2` (the polarity fixed point, the qubit dimension's unsigned magnitude)
- maxval at `1/4 = (1/2)²` (the cardioid cusp, the Theorem-1 ceiling, the operator-space mirror partner of the discriminant 4 = 2²)

Coda from `reflections/ON_THE_HALF.md:39-45` (Tom 2026-05-07):

> "1/4 is the half's quadratic shadow … two readings, one parabola; the half is the axis, the quarter is the height."

This is the meta-pattern that the inheritance tree formalises: every Tier1Derived anchor in Pi2KB either IS one of these two values or factors algebraically through them.

---

## Cross-references

- Foundational claims: `compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs`
- Wiring centre: `compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBase.cs`
- F1 root: `compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs`
- F1 proof: `docs/proofs/MIRROR_SYMMETRY_PROOF.md`
- Theorem-1 proof: `docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md`
- Coda & synthesis: `reflections/ON_THE_HALF.md`, `reflections/ON_BOTH_SIDES_OF_THE_MIRROR.md`
- Bridge framing: `docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md`
- F-formula registry: `docs/ANALYTICAL_FORMULAS.md`
- Symmetry-axis inventory (parallel structure on the dynamics-symmetry side): `docs/SYMMETRY_FAMILY_INVENTORY.md`

## Z2Axis classification (added 2026-05-24, updated 2026-05-28)

Every Pi²-Inheritance Claim carries a `Z2Axis` enum classification. Counts are live
via `rcpsi inspect PolarityCubeMap`; the bit_a side has grown from 1 to 7 since the
original 2026-05-24 snapshot (Welle 7/15).

| Z2Axis value | Meaning | Count |
|--------------|---------|-------|
| BitB | Π²_Z = X⊗N axis (F1² family, n_Y + n_Z parity) | 67 |
| BitA | Π²_X = Z⊗N axis (F61 family, n_X + n_Y parity) | 7 (F61, F38BitA, F39BitA, F63BitARef, ZMirrorBitA, F108Part2, F112-X) |
| Klein2 | Uses both Π²_Z and Π²_X axes (Klein-Vierergruppe) | KleinFourCellClaim + F79 / F89F88aKlein / Pi2KleinBilinearTable |
| YParity | Term-level Y-parity refinement (k≥3) | F102 + the F87 Z₂³ family |
| Cubic3 | Full Z₂³ (Klein-eight) | KleinEightCellClaim + ZGlobalMirrorRefinement |
| NotApplicable | Foundation Claims (Pi2InvolutionClaim, etc.) | live |

### bit_a-twin coverage and its completion

Of the 67 BitB claims, the bit_a-twin slots are: **7 Filled** (typed bit_a sibling
wired), **8 BitBSpecific** (no twin by construction: amplitude damping F1T1 / F82 /
F84, plus F91, F93, F108Part3-Y, F113, F112-Y), **17 CoveredByHadamardDuality** (the
Absorption-Theorem descendant family F33 / F50 / F55 / F64-F68 / F74, plus F83 the
anti-fraction, the operator-space Π-decomposition family F49 / F49b / F80 / F81, and the
Lindblad-spectral Absorption family F3 / F43 / F44, now typed: the bit_a twin is a
corollary of the duality below, no bespoke twin owed), and **35 NeedsDerivation** (open).

These open slots are not independent problems. They are governed by one proven principle,
the **global Hadamard X↔Z duality** (the Klein-V₄ outer automorphism,
[`PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md`](proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md):
`Q_zx · Π_Z · Q_zx⁻¹ = Π_X`), consolidated for the twin classification in
[`PROOF_BIT_A_TWIN_VIA_HADAMARD.md`](proofs/PROOF_BIT_A_TWIN_VIA_HADAMARD.md). The
enabler is the **Absorption Theorem's H-independence** (proven for any Hermitian H,
2026-05-28): since only the dissipator sets Re(λ), swapping the dissipator axis Z↔X
by the global Hadamard transports a bit_b result to its bit_a twin.

- **Collapsible by the duality** (the majority, status `CoveredByHadamardDuality`):
  any claim reducing to Π/L spectrum, eigenspace, or operator-identity content. The
  nine BitB members of the Absorption-descendant family (F33, F50, F55, F64-F68, F74),
  plus F83 (the anti-fraction, an operator-space Frobenius-norm statement on the 4^N
  Pauli space), plus F49 / F49b / F80 / F81 (the operator-space Π-decomposition family:
  residual and dissipator Frobenius norms, Spec(M), and the Π-conjugation identity), plus
  F3 / F43 / F44 (the Lindblad-spectral Absorption family: decay-rate bounds, sector SFF
  pairing, and the palindromic eigenvalue-pair rate identity), are
  now typed with this status; F89 is itself an AT-descendant but sits on the Klein2
  axis, so it has no bit_a-twin slot. The operator-identity and spectral/mirror claims
  are corollaries of the same theorem and form the next batch to type.
- **Bespoke-operator residue** (a handful): claims built on bespoke operators
  (Π_5bilinear, F108-style) where the operator-space Hadamard does not transport
  directly and a deeper Hilbert-space Hadamard is needed. These keep
  `NeedsDerivation`.
- **Lift caveat:** only {I, Q_zx (Hadamard)} lifts to a Hilbert-space unitary, so only
  Z↔X transports Lindblad form; D and Q_yx are operator-space-only.

The `PolarityCubeMap` typed Claim aggregates this inventory and exposes it via
`rcpsi inspect PolarityCubeMap`.

Source: `compute/RCPsiSquared.Core/Symmetry/Z2Axis.cs`,
`IZ2AxisClaim.cs`, `PolarityCubeMap.cs`.
