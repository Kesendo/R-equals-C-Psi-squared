using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>F87's trichotomy {truly, soft, hard} inherits from the Pi2-Foundation
/// at the number level via F1's residual M.
///
/// <para><b>F1 residual M as discriminator</b>: F87 classifies via the F1 palindrome
/// residual <c>M = Π·L·Π⁻¹ + L + 2σ·I</c>. The "2" coefficient in F1's M is already
/// typed via <see cref="F1Pi2Inheritance"/> as <c>a_0</c> on the Pi2 dyadic ladder.
/// F87 inherits this anchor transitively.</para>
///
/// <para><b>F87 ≠ Klein cells</b>: F87 (trichotomy 14/19/3 from Pauli-pair
/// combinatorics, commit 81caf67) was developed before the
/// <see cref="KleinFourCellClaim"/> 4-cell formalism. The two structures are distinct,
/// even though both involve Π². Earlier framings of F87 as a "4-way Π²-class
/// refinement aligned with KleinFour" were a post-hoc overlay; this claim no longer
/// asserts that edge. Klein cells stay typed via their own primitives.</para>
///
/// <para>Tier1Derived: pure inheritance composition. F87 itself is Tier1Derived
/// (<see cref="F87TrichotomyClassification"/>); F1's M residual is Tier1Derived
/// in <c>docs/proofs/MIRROR_SYMMETRY_PROOF.md</c>.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F87 +
/// <c>experiments/V_EFFECT_FINE_STRUCTURE.md</c> +
/// <c>compute/RCPsiSquared.Diagnostics/F87/F87TrichotomyClassification.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs</c>.</para></summary>
public sealed class F87Pi2Inheritance : Claim
{
    private readonly F1Pi2Inheritance _f1Inheritance;

    /// <summary>The "2" coefficient transitively inherited from F1's residual via
    /// <see cref="F1Pi2Inheritance.TwoFactor"/>. F87's discriminator <c>M</c> uses this
    /// "2" inside the F1 closed form; F87 inherits its number-level constant from
    /// the Pi2 dyadic ladder via F1.</summary>
    public double TransitivelyInheritedTwoFactor => _f1Inheritance.TwoFactor;

    public F87Pi2Inheritance(F1Pi2Inheritance f1Inheritance)
        : base("F87 trichotomy inherits from Pi2-Foundation: discriminator via F1 (2 = a_0)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F87 + " +
               "experiments/V_EFFECT_FINE_STRUCTURE.md + " +
               "compute/RCPsiSquared.Diagnostics/F87/F87TrichotomyClassification.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs")
    {
        _f1Inheritance = f1Inheritance ?? throw new ArgumentNullException(nameof(f1Inheritance));
    }

    public override string DisplayName =>
        "F87 ←(F1 residual)← Pi2-Foundation";

    public override string Summary =>
        $"F87 trichotomy {{truly, soft, hard}} discriminated via F1 residual M (2 = a_0 = {TransitivelyInheritedTwoFactor}) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F87 closed form",
                summary: "trichotomy via F1 residual: truly iff ‖M‖ < ε; soft iff M ≠ 0 but spectrum pairs; hard iff no pairing");
            yield return new InspectableNode("number-level inheritance via F1",
                summary: $"F87's discriminator M uses F1's 2σ·I shift; the 2 = a_0 = {TransitivelyInheritedTwoFactor} comes from Pi2 dyadic ladder");
            yield return InspectableNode.RealScalar("TransitivelyInheritedTwoFactor (= a_0 via F1)", TransitivelyInheritedTwoFactor);
            yield return new InspectableNode("F87 hardware confirmation",
                summary: "Marrakesh d7mjnjjaq2pc73a1pk4g (2026-04-26) Δ(soft − truly) = −0.722; Kingston regime-uniformity (2026-05-05)");
            yield return new InspectableNode("V-Effect counts",
                summary: "36-enum N=3 → 14/19/3 (truly/soft/hard); 120-enum N=3..5 → 15/46/59 N-stable (combinatorial proof commit 81caf67)");
            yield return new InspectableNode("F87 vs Klein cells",
                summary: "distinct structures: F87 trichotomy (14/19/3, 36-enum) developed pre-Klein from Pauli-pair combinatorics; KleinFour 4-cell (Pp/Pm/Mp/Mm) is a separate Π²-eigenspace decomposition; no typed inheritance edge between them");
        }
    }
}
