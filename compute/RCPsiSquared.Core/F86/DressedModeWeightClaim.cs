using RCPsiSquared.Core.Inspection;

using RCPsiSquared.Core.Knowledge;
namespace RCPsiSquared.Core.F86;

/// <summary>F86 Tier-1 derived structural fingerprint: at Q = Q_peak the Dicke probe sits
/// approximately 99 % in dressed (H-mixed, complex-eigenvalue) modes; at Q = 20 (plateau,
/// far past Q_peak) it sits ≈ 31 % in dressed modes.
///
/// <para>This is what makes Q_peak a generalised exceptional-point resonance condition: the
/// J-derivative ∂S/∂J peaks where the probe weight has been pulled off the pure-rate ladder
/// onto the first complex-conjugate eigenvalue pair just past the EP. The 99 % at peak vs
/// 31 % at plateau is the structural signature documented in <c>ANALYTICAL_FORMULAS.md</c>
/// F86 EP-based structural derivation paragraph.</para>
/// </summary>
public sealed class DressedModeWeightClaim : Claim
{
    public double WeightAtQPeak { get; } = 0.99;
    public double WeightAtPlateau { get; } = 0.31;
    public double PlateauQ { get; } = 20.0;

    public DressedModeWeightClaim()
        : base("dressed-mode probe weight at Q_peak",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F86 EP-based structural derivation")
    { }

    public override string DisplayName => "Dicke probe → dressed modes: 99% at Q_peak, 31% at Q=20";

    public override string Summary =>
        $"at Q_peak: {WeightAtQPeak:P0}; at Q={PlateauQ:F0} plateau: {WeightAtPlateau:P0} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("dressed-mode weight at Q_peak", WeightAtQPeak, "P1");
            yield return InspectableNode.RealScalar($"dressed-mode weight at Q={PlateauQ}", WeightAtPlateau, "P1");
            yield return new InspectableNode("interpretation",
                summary: "probe pulled off pure-rate ladder onto complex-conjugate pair just past EP — Q_peak is a generalised EP resonance");
        }
    }
}
