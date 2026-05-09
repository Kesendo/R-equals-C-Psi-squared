using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F74ChromaticityPi2Inheritance"/>:
/// F74 chromaticity c(n, N) = min(n, N−1−n) + 1. One typed parent edge:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c>
///         (RateCoefficient in 2γ₀·HD pure rates).</item>
/// </list>
///
/// <para>Tier consistency: F74 is Tier 1 combinatorial; verified N=3..7.</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>.</para></summary>
public static class F74ChromaticityPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF74ChromaticityPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F74ChromaticityPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            return new F74ChromaticityPi2Inheritance(ladder);
        });
}
