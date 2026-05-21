using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="LindbladAbsorptionMatchAtSixtyDegreesClaim"/>: at
/// Q = √3 the Lindblad 2×2 sub-block eigenvalue magnitude |λ_±| = γ₀·√(1+Q²) equals the
/// Absorption Theorem single-site rate 2γ₀, and the F95 angle θ = arctan(Q) lands on the
/// canonical Niven angle 60°. Three parent edges:
///
/// <list type="bullet">
///   <item><see cref="F95AngleAtQuadraticZeroPi2Inheritance"/>: the Lindblad 2×2 angle law
///         θ = arctan(Q).</item>
///   <item><see cref="AbsorptionTheoremClaim"/>: the single-site rate α = 2γ₀.</item>
///   <item><see cref="CanonicalTrigAnchorPi2Inheritance"/>: the canonical Niven angles, of
///         which 60° is the match point.</item>
/// </list>
///
/// <para>The registry instance is built from the registered parents. The standalone
/// <see cref="LindbladAbsorptionMatchAtSixtyDegreesClaim.Build"/> and <c>Shared</c> factories
/// (used by F86KnowledgeBase) construct their own parent chain and are independent of the
/// registry.</para>
///
/// <para>Tier consistency: all four Tier1Derived.</para>
///
/// <para>Requires upstream registrations:
/// <see cref="F95AngleAtQuadraticZeroPi2InheritanceRegistration.RegisterF95AngleAtQuadraticZeroPi2Inheritance"/> +
/// <see cref="AbsorptionTheoremClaimRegistration.RegisterAbsorptionTheoremClaim"/> +
/// <see cref="CanonicalTrigAnchorPi2InheritanceRegistration.RegisterCanonicalTrigAnchorPi2Inheritance"/>.</para></summary>
public static class LindbladAbsorptionMatchAtSixtyDegreesClaimRegistration
{
    public static ClaimRegistryBuilder RegisterLindbladAbsorptionMatchAtSixtyDegreesClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<LindbladAbsorptionMatchAtSixtyDegreesClaim>(b =>
        {
            var f95 = b.Get<F95AngleAtQuadraticZeroPi2Inheritance>();
            var absorption = b.Get<AbsorptionTheoremClaim>();
            var canonicalTrig = b.Get<CanonicalTrigAnchorPi2Inheritance>();
            return new LindbladAbsorptionMatchAtSixtyDegreesClaim(f95, absorption, canonicalTrig);
        });
}
