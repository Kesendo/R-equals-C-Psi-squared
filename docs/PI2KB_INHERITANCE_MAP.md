# Pi2 Knowledge Base: Inheritance Map

**Question Tom 2026-05-16:** *"Ich denke das durch Vererbung diese Anker vielleicht auf irgendeine Art zusammenhängen, irgendwo müssen sich die Anker vererben."*

**Answer:** yes, the structural anchors of the R=CΨ² framework are organised as a lineage tree rooted in `compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs`. This document is the visual map of that lineage, the honest split between its typed and its semantic edges, and the live pointers into the Object Manager.

**Date:** 2026-05-16 (created), last refreshed 2026-07-16 (the change history lives in git)

---

## The tree (the semantic lineage)

```
PolynomialFoundationClaim                  Pi2KnowledgeBaseClaims.cs
│   d²−2d=0 ↔ R=CΨ² ; minimum-memory equation
│   The trunk; all four branches root here.
│
├── number-branch (1/d = 1/2 lineage)
│   │
│   QubitDimensionalAnchorClaim
│   │   1/d = 1/2 at d=2 (the unsigned magnitude basis)
│   │
│   ├── BilinearApexClaim
│   │   │  1/2 = argmax of p·(1−p)
│   │   │
│   │   ├── QuarterAsBilinearMaxvalClaim
│   │   │      1/4 = (1/2)² = maxval of p·(1−p)
│   │   │
│   │   └── ArgmaxMaxvalPairClaim
│   │          closure: (1/2 argmax, 1/4 maxval) inseparable invariant
│   │
│   ├── HalfAsStructuralFixedPointClaim
│   │      three faces of 1/2 close (bridge, horizon, substrate)
│   │
│   └── Pi2DyadicLadderClaim               Pi2DyadicLadderClaim.cs
│       │   a_n = 2^(1−n) ladder rungs
│       │
│       ├── PolynomialDiscriminantAnchorClaim   [typed]
│       │      a_{−1} = 4 = polynomial discriminant of d²−2d=0
│       │      mirror partner: a_3 = 1/4 (the QuarterAsBilinearMaxval); closure 4·(1/4)=1
│       │
│       ├── Pi2OperatorSpaceMirrorClaim
│       │
│       ├── AbsorptionTheoremClaim              [typed]
│       │      α = 2γ₀·⟨n_XY⟩  (quantum decay rate)
│       │      └── UniversalCarrierClaim        [typed]
│       │
│       └── the F-leaves                        [typed: most take the ladder; see below]
│              60 F-numbered *Pi2Inheritance.cs files, 56 distinct F tokens
│              from F1 to F97 (counting F2b/F49b/F49c as their own tokens;
│              the file listing IS the exact inventory, the numbering has gaps),
│              plus 4 non-F leaves (CanonicalTrigAnchor, DickeSuperpositionQuarter,
│              QubitNecessity, XGlobalEigenstateMirror). The bit_a-axis leaves
│              (F38BitAInvolutionInheritance, F39DetPiBitAInheritance,
│              F63BitAReference, ZGlobalEigenstateMirrorBitAInheritance) and
│              the mixed-axis F108 and F112/F113 families (per-part axes;
│              F113 itself is bit_b) live outside the *Pi2Inheritance naming
│              pattern and are counted in the PolarityCubeMap inventory below.
│
├── angle-branch (90° lineage)
│   │
│   NinetyDegreeMirrorMemoryClaim
│   │   90° rotation back to the mirror (the i in F80's 2i)
│   │
│   ├── Pi2InvolutionClaim                      [typed: ctor injects F1PalindromeIdentity]
│   │   │   Π² commutes with L (squared from F1)
│   │   │
│   │   └── KleinFourCellClaim
│   │          F88a two-axis (Π²_Z, Π²_X) → 4 Klein cells
│   │
│   └── Pi2I4MemoryLoopClaim                    Pi2I4MemoryLoopClaim.cs
│          i⁴ = 1 (Z₄ closure of 90° rotations)
│
├── polarity branch
│   PolarityLayerOriginClaim
│       +0 / 0 / −0 polarity layer at d=0 → ±0.5 pair at d=2
│
└── perspectival branch
    TwoReadingsClaim                            TwoReadingsClaim.cs
        any object on d²=4^N admits two coordinate readings of one
        underlying object. Names seven layers of this pattern:
        number/angle of d, argmax/maxval of p·(1−p), M/Π·M·Π⁻¹ (F81),
        bra/ket of ρ, inside/outside Q = J/γ₀, classical/quantum Lese-Modus,
        inter-sectoral wave ("we are the standing wave between").
```

The number and angle branches **reconverge** at the F-leaves (most cite the Dyadic ladder, many also the I4 Memory Loop; a handful anchor elsewhere, see below); the ladder's rungs hold the 1/2-lineage values.

---

## How to read the tree: semantic lineage vs typed edges

The tree above is the SEMANTIC map (which anchor generates which). Only a
subset of its edges are typed ctor injections, marked `[typed]`: a child
Claim's constructor accepts the parent Claim and stores it, and
`Pi2KnowledgeBase.cs` does the wiring centrally. The unmarked edges are
documented lineage (docstrings, anchor strings, proofs), deliberate but not
enforced by the type system.

The typed-graph ROOTS, by the wiring's own definition (claims constructed
with no Claim-typed ctor parents in `Pi2KnowledgeBase.cs`), are twelve:
PolynomialFoundation, QubitDimensionalAnchor, F1PalindromeIdentity,
KleinFourCell, BilinearApex, QuarterAsBilinearMaxval, ArgmaxMaxvalPair,
HalfIntegerMirror (takes only `int N`), HalfAsStructuralFixedPoint,
NinetyDegreeMirrorMemory, PolarityLayerOrigin, Pi2DyadicLadder. The claims
built WITH typed parents there are Pi2Involution (from F1),
PolynomialDiscriminantAnchor, AbsorptionTheorem, UniversalCarrier, and
CanonicalTrigAnchor, plus the local intermediates
DickeSuperpositionQuarter and the F98 long-time claim feeding
CanonicalTrigAnchor. Most F-leaves take `Pi2DyadicLadderClaim` in their
ctor (and often I4 Memory Loop, Half, ...), wired in
`compute/RCPsiSquared.Runtime/PolarityArchitecture/`; seven anchor
elsewhere: F91/F92/F93 on I4 Memory Loop alone, F95/F97 on
Half/Quarter/NinetyDegree/Foundation, F79 on
KleinFourCell/OperatorSpaceMirror/F1, F86F71 on upstream F-claims.

The SEMANTIC trunk remains three-rooted: `PolynomialFoundationClaim` (the
trunk of all four branches), `F1PalindromeIdentity`
(`compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs`, the Π branch via
Π²·L·Π⁻² = L), and `ChiralKClaim`
(`compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs`, the chiral/sublattice
root K = diag((−1)^ℓ), spectrum inversion at the Hamiltonian level). None
derives from another; ChiralK is wired in the wider claim registry, not in
`Pi2KnowledgeBase`.

---

## Open edges (current)

1. **F71 → HalfIntegerMirror**: prose by design. `F71MirrorSymmetryPi2Inheritance`
   is N-universal; `HalfIntegerMirrorClaim` takes `int N` at construction, so
   the edge cannot be a ctor parent without parameterising F71 by N or a
   per-N adapter. The ctor takes only `Pi2DyadicLadderClaim`.
2. **F77 → F75 (Taylor-limit mother)**: held as registration discard. The
   F75-chain construction (`F66 → F65 → F75`) is multi-step; F77's ctor
   takes the ladder and `HalfAsStructuralFixedPointClaim` only.
3. **The F86 1/4 anchors without graph edges**: `BlockCpsiTrajectory`
   hard-codes the crossing threshold `0.25` as a method default and cites
   the roadmap only in its anchor string; `BlockCoherenceContent` is a
   static class (no ctor possible) holding `Quarter = 0.25` as a const with
   `<see cref>` pointers. Neither carries a typed 1/4 parent; the typed F86
   legs that exist are `IbmBlockCpsiHardwareTable` and
   `F86TPeakPi2Inheritance` (whose ctor takes `QuarterAsBilinearMaxvalClaim`).

---

## The anchor-side sweep, realized

The complement to top-down parent-finding ("start from each anchor, find
where its constant is consumed without a typed edge") exists as live
objects: `F99AnchorMap` aggregates the `IF99AnchorBearing` claims by anchor
value {0, 1/8, 1/4, 3/8, 1/2} with explicit gap detection (`GapAnchors`),
`DickeAnchorMap` does the same for the Dicke anchors, and
`RegistryWiringAuditTests` (Runtime.Tests/AnchorAudit) audits the wiring.
The dominant framing since is the cubic classification below: an un-typed
anchor consumption shows up as an open twin slot.

---

## The single numerical structure underneath

All anchors collapse to **one quadratic parabola** `p · (1 − p)`:

- argmax at `p = 1/2` (the polarity fixed point, the qubit dimension's unsigned magnitude)
- maxval at `1/4 = (1/2)²` (the cardioid cusp, the Theorem-1 ceiling, the operator-space mirror partner of the discriminant 4 = 2²)

Coda from `reflections/ON_THE_HALF.md` (Tom 2026-05-07):

> "1/4 is the half's quadratic shadow … two readings, one parabola; the half is the axis, the quarter is the height."

This is the meta-pattern the inheritance tree formalises: every Tier1Derived anchor in Pi2KB either IS one of these two values or factors algebraically through them.

---

## Z2Axis classification

Every Pi²-Inheritance Claim carries a `Z2Axis` enum classification. The
live inventory is `PolarityCubeMap`
(`compute/RCPsiSquared.Core/Symmetry/PolarityCubeMap.cs`, registered via
`PolarityCubeMapRegistration` into the world registry): inspect it with

    dotnet run --project compute/RCPsiSquared.Cli -- inspect --claim PolarityCubeMap

(it also appears under `--root world`; the `--root pi2` view shows only the
4-claim local core). The inventory spans the Core + Runtime registry;
Diagnostics-hosted leaves (e.g. F87Pi2Inheritance) are excluded by design.
Live counts (2026-07-16, 89 claims in the inventory):

| Z2Axis value | Meaning | Count |
|--------------|---------|-------|
| BitB | Π²_Z = X⊗N axis (F1² family, n_Y + n_Z parity) | 67 |
| BitA | Π²_X = Z⊗N axis (F61 family, n_X + n_Y parity) | 7 |
| Klein2 | Uses both Π²_Z and Π²_X axes (Klein-Vierergruppe) | 5 |
| YParity | Term-level Y-parity refinement (k≥3) | 8 |
| Cubic3 | Full Z₂³ (Klein-eight) | 1 |
| NotApplicable | Foundation-layer claim without an axis | 1 |

### bit_a-twin coverage

Of the 67 BitB claims (twin coverage 10.4%): **7 Filled** (typed bit_a
sibling wired), **8 BitBSpecific** (no twin by construction: amplitude
damping F1T1 / F82 / F84, plus F91, F93, F108Part3-Y, F113, F112-Y),
**19 CoveredByHadamardDuality**, and **33 NeedsDerivation** (the open twin
slots).

The open slots are not independent problems. They are governed by one
proven principle, the **global Hadamard X↔Z duality** (the Klein-V₄ outer
automorphism,
[`PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md`](proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md):
`Q_zx · Π_Z · Q_zx⁻¹ = Π_X`), consolidated for the twin classification in
[`PROOF_BIT_A_TWIN_VIA_HADAMARD.md`](proofs/PROOF_BIT_A_TWIN_VIA_HADAMARD.md).
The enabler is the **Absorption Theorem's H-independence** (proven for any
Hermitian H): since only the dissipator sets Re(λ), swapping the dissipator
axis Z↔X by the global Hadamard transports a bit_b result to its bit_a twin.

- **Collapsible by the duality** (`CoveredByHadamardDuality`): any claim
  reducing to Π/L spectrum, eigenspace, or operator-identity content: the
  Absorption-descendant family (F33, F50, F55, F64-F68, F74), F83, the
  operator-space Π-decomposition family (F49 / F49b / F80 / F81), the
  Lindblad-spectral family (F3 / F43 / F44), F49c, and the k-body roll-up
  F85. F89 is itself an AT-descendant but sits on the Klein2 axis, so it
  has no bit_a-twin slot.
- **Bespoke-operator residue** (`NeedsDerivation`): claims built on bespoke
  operators (Π_5bilinear, F108-style) where the operator-space Hadamard
  does not transport directly and a deeper Hilbert-space Hadamard is needed.
- **Lift caveat:** only {I, Q_zx (Hadamard)} lifts to a Hilbert-space
  unitary, so only Z↔X transports Lindblad form; D and Q_yx are
  operator-space-only.

Source: `compute/RCPsiSquared.Core/Symmetry/Z2Axis.cs`, `IZ2AxisClaim.cs`,
`PolarityCubeMap.cs`.

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
