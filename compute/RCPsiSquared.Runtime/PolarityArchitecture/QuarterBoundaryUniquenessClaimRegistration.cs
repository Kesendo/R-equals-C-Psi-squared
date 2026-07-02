using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="QuarterBoundaryUniquenessClaim"/>. One parent edge:
/// <see cref="QuarterAsBilinearMaxvalClaim"/> (the ¼ value the uniqueness argument singles out),
/// registered earlier in <c>BuildDefault</c> via the Pi2 family.
///
/// <para>Requires: RegisterQuarterAsBilinearMaxvalClaim (Pi2 family).</para></summary>
public static class QuarterBoundaryUniquenessClaimRegistration
{
    public static ClaimRegistryBuilder RegisterQuarterBoundaryUniquenessClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<QuarterBoundaryUniquenessClaim>(b =>
            new QuarterBoundaryUniquenessClaim(b.Get<QuarterAsBilinearMaxvalClaim>()));
}
