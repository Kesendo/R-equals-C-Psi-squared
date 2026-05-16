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
│              F81, F82, F83, F84, F85, F86 family, F87, F88, F89, F90,
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
│   │          F88 two-axis (Π²_Z, Π²_X) → 4 Klein cells
│   │
│   └── Pi2I4MemoryLoopClaim                Pi2I4MemoryLoopClaim.cs
│          i⁴ = 1 (Z₄ closure of 90° rotations)
│
└── polarity branch
    PolarityLayerOriginClaim                :516
        +0 / 0 / −0 polarity layer at d=0 → ±0.5 pair at d=2
```

Both branches **reconverge** at `Pi2DyadicLadderClaim` (rungs hold 1/2 lineage values) and at the F-leaves (each F-claim cites both Dyadic + I4 Memory Loop).

---

## Tier-1 roots

Two roots in the typed graph (no inbound edges from other Claims):

1. **`PolynomialFoundationClaim`** (`Pi2KnowledgeBaseClaims.cs:441`) — semantic root of the number-anchor lineage. Its docstring states: *"This is the trunk that generates both framework anchors: `HalfAsStructuralFixedPointClaim` (1/2 number-anchor) and `NinetyDegreeMirrorMemoryClaim` (90° angle-anchor)."* Only markdown citations.

2. **`F1PalindromeIdentity`** (`compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs:24`) — root of the angle/Π-branch via its squaring step `Π² · L · Π⁻² = L`. Now explicitly injected into `Pi2InvolutionClaim` as a typed ctor parent (2026-05-16); before that the edge was prose only.

---

## Remaining prose-only edges (asserted but not yet typed)

Each line is an edge that EXISTS in narrative form (docstring, proof, reflection) but is NOT yet a typed Claim-to-Claim ctor injection.

| Edge | Asserted at | Why not typed (yet) |
|---|---|---|
| **F71 spatial mirror → F1 / Π cluster** | `MIRROR_SYMMETRY_PROOF.md`, `PROOF_C1_MIRROR_SYMMETRY.md` | F71MirrorSymmetryPi2Inheritance attaches to `Pi2DyadicLadderClaim` directly; no F1-sibling edge |
| **K chiral / sublattice (BDI class)** | `compute/RCPsiSquared.Core/Symmetry/ChiralK.cs`, memory `project_chiral_partnership` | No Claim wrapper for K; only a helper class |
| **F86, F57, F64 as siblings of 1/4** | `docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md`, `docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md` | Each F-claim attaches independently to the ladder; no typed "1/4-sibling" edge |
| **"Two readings" anchor** | `reflections/ON_BOTH_SIDES_OF_THE_MIRROR.md`, `reflections/ON_THE_HALF.md` | No Claim wrapper for the perspectival reading itself |
| **F1 → F71 spatial → F81 chain** | `docs/proofs/PROOF_F81_PI_CONJUGATION_OF_M.md` | Each leaf attaches directly to `Pi2DyadicLadderClaim`; F-to-F lineage is narrative only |

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
