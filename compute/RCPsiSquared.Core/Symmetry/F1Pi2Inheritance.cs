using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F1's palindrome closed form <c>Π · L · Π⁻¹ = −L − 2Σγ · I</c> has its
/// "<c>2</c>" coefficient sitting on the Pi2-Foundation: it is exactly <c>a_0</c> on
/// the dyadic halving ladder, the qubit dimension d. The shift "2Σγ·I" is therefore
/// "<c>d · Σγ · I</c>" — the Σγ accumulates the per-site dephasing rates, and the
/// "d" counts how the F1 derivation stacks two Π applications.
///
/// <para>The "+L" and "−L" form a Z₂ inversion in the Liouvillian operator algebra;
/// applying Π once gives <c>−L − 2σ·I</c>, applying it twice closes back to <c>+L</c>.
/// That two-step closure is the F1 palindrome reading: Π is involutive on the
/// algebra modulo the constant shift <c>2σ·I</c>, where the "2" is the qubit
/// dimension.</para>
///
/// <para>This is not a new identity — F1 is Tier1Derived in
/// <c>docs/proofs/MIRROR_SYMMETRY_PROOF.md</c> and verified at N=8 (54,118 eigenvalues,
/// zero exceptions). What this claim makes typed is the <b>inheritance</b>: the
/// constant 2 in F1's closed form is not a free parameter; it is <c>a_0 = d</c> on
/// the Pi2 dyadic ladder.</para>
///
/// <para>Tier1Derived: pure composition. Anchors:
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F1 +
/// <c>docs/proofs/MIRROR_SYMMETRY_PROOF.md</c> +
/// <c>compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F1Pi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;

    /// <summary>The "2" coefficient in F1's "−L − 2Σγ·I" closed form. Exactly equal to
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(0) = <c>a_0</c> = d (qubit dimension).</summary>
    public double TwoFactor => _ladder.Term(0);

    public F1Pi2Inheritance(Pi2DyadicLadderClaim ladder)
        : base("F1 palindrome's 2 coefficient inherits from Pi2-Foundation (2 = a_0)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F1 + docs/proofs/MIRROR_SYMMETRY_PROOF.md + " +
               "compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
    }

    public override string DisplayName =>
        "F1's 2 coefficient as Pi2-Foundation a_0";

    public override string Summary =>
        $"Π·L·Π⁻¹ = −L − 2Σγ·I: the 2 coefficient = a_0 = {TwoFactor} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F1 closed form",
                summary: "Π·L·Π⁻¹ = −L − 2Σγ·I (Tier1Derived in MIRROR_SYMMETRY_PROOF)");
            yield return InspectableNode.RealScalar("TwoFactor (= a_0 = d)", TwoFactor);
            yield return new InspectableNode("Z₂ inversion",
                summary: "+L and −L − 2σ·I are Π-partners; two Π applications close back; the 2 is the qubit dimension d, not a free parameter");
        }
    }
}
