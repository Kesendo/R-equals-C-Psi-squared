using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="CanonicalTrigAnchorPi2Inheritance"/> (F99): the
/// F86b α(γ) = (1−γ²)/2 closed form evaluated at the canonical trig angles
/// {0°, 30°, 45°, 60°, 90°}, producing the five Pi2 dyadic anchors {0, 1/8, 1/4, 3/8, 1/2}.
/// Three parent edges:
///
/// <list type="bullet">
///   <item><see cref="HalfAsStructuralFixedPointClaim"/>: the 90° anchor (Generic Dicke).</item>
///   <item><see cref="QuarterAsBilinearMaxvalClaim"/>: the 45° anchor (non-uniform Dicke).</item>
///   <item><see cref="KIntermediateAsymptoteQuarterInheritance"/>: the F98 long-time
///         companion, the dynamic 3/8 to 1/4 bridge.</item>
/// </list>
///
/// <para>Tier consistency: all four Tier1Derived.</para>
///
/// <para>Requires upstream registrations: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (HalfAsStructuralFixedPointClaim, QuarterAsBilinearMaxvalClaim) +
/// <see cref="KIntermediateAsymptoteQuarterInheritanceRegistration.RegisterKIntermediateAsymptoteQuarterInheritance"/>.</para></summary>
public static class CanonicalTrigAnchorPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterCanonicalTrigAnchorPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<CanonicalTrigAnchorPi2Inheritance>(b =>
        {
            var half = b.Get<HalfAsStructuralFixedPointClaim>();
            var quarter = b.Get<QuarterAsBilinearMaxvalClaim>();
            var f98 = b.Get<KIntermediateAsymptoteQuarterInheritance>();
            return new CanonicalTrigAnchorPi2Inheritance(half, quarter, f98);
        });
}
