using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 Statement 1 partial: <c>Q_EP = 2/g_eff</c> — exceptional-point Q at which
/// the slowest 2-level eigenvalue pair coalesces.
///
/// <para>Tier 1 derived in the 2-level reduction; the value matches the full block-L EP
/// position bit-exactly when g_eff is taken as the SVD-top inter-channel singular value σ_0.
/// For c=2, σ_0 → 2√2 asymptotically (numerical witness across N=5..8).</para>
/// </summary>
public sealed class QEpLaw : F86Claim
{
    public double GEff { get; }
    public double Value { get; }

    public QEpLaw(double gEff)
        : base("Q_EP law", Tier.Tier1Derived, "docs/ANALYTICAL_FORMULAS.md F86 Statement 1")
    {
        GEff = gEff;
        Value = EpAlgebra.QEp(gEff);
    }

    public override string DisplayName => "Q_EP = 2/g_eff";
    public override string Summary => $"= {Value:G6} (g_eff = {GEff:G6}, {Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("g_eff", GEff, "G6");
            yield return InspectableNode.RealScalar("Q_EP", Value, "G6");
        }
    }

    public override InspectablePayload Payload =>
        new InspectablePayload.Real("Q_EP", Value, "G6");
}
