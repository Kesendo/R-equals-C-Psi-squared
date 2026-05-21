using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="KIntermediateAsymptoteQuarterInheritance"/> (F98):
/// the KIntermediate Dicke long-time Π²-odd asymptote (N+2)/[4(N+1)] bridging the F86b 3/8
/// static anchor to the 1/4 universal boundary. Three parent edges:
///
/// <list type="bullet">
///   <item><see cref="QuarterAsBilinearMaxvalClaim"/>: the 1/4 asymptote target.</item>
///   <item><see cref="HalfAsStructuralFixedPointClaim"/>: the 1/4 = (1/2)² polarity-squared
///         algebra.</item>
///   <item><see cref="DickeSuperpositionQuarterPi2Inheritance"/>: the static-side companion
///         (the C_block = 1/4 ceiling for the same Dicke superposition family).</item>
/// </list>
///
/// <para>Tier consistency: all four Tier1Derived.</para>
///
/// <para>Requires upstream registrations: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (QuarterAsBilinearMaxvalClaim, HalfAsStructuralFixedPointClaim) +
/// <see cref="DickeSuperpositionQuarterPi2InheritanceRegistration.RegisterDickeSuperpositionQuarterPi2Inheritance"/>.</para></summary>
public static class KIntermediateAsymptoteQuarterInheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterKIntermediateAsymptoteQuarterInheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<KIntermediateAsymptoteQuarterInheritance>(b =>
        {
            var quarter = b.Get<QuarterAsBilinearMaxvalClaim>();
            var half = b.Get<HalfAsStructuralFixedPointClaim>();
            var staticSide = b.Get<DickeSuperpositionQuarterPi2Inheritance>();
            return new KIntermediateAsymptoteQuarterInheritance(quarter, half, staticSide);
        });
}
