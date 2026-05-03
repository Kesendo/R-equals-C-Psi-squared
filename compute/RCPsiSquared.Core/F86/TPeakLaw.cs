using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 Statement 1: <c>t_peak = 1/(4γ₀)</c> — universal e-folding time at the
/// slowest channel pair, independent of c, N, n, bond position. Bit-exact verified against
/// full block-L numerics across all tested cases.
///
/// <para>This is THE most-derived F86 claim — the analytic 2-level reduction's answer that
/// matches the full Liouvillian numerics exactly. Construct from γ₀ and the value falls out
/// trivially via <see cref="EpAlgebra.TPeak"/>.</para>
/// </summary>
public sealed class TPeakLaw : F86Claim
{
    public double GammaZero { get; }
    public double Value { get; }

    public TPeakLaw(double gammaZero)
        : base("t_peak law", Tier.Tier1Derived, "docs/ANALYTICAL_FORMULAS.md F86 Statement 1")
    {
        GammaZero = gammaZero;
        Value = EpAlgebra.TPeak(gammaZero);
    }

    public override string DisplayName => "t_peak = 1/(4γ₀)";
    public override string Summary => $"= {Value:G6} ({Tier.Label()}, universal across c, N, n, bond)";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("γ₀", GammaZero, "G6");
            yield return InspectableNode.RealScalar("t_peak", Value, "G6");
        }
    }

    public override InspectablePayload Payload =>
        new InspectablePayload.Real("t_peak", Value, "G6");
}
