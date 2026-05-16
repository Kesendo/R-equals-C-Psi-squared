using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F25CPsiBellPlusPi2Inheritance"/>:
/// F25's <c>CΨ(t) = f(1+f²)/6</c> Bell+ Z-dephasing closed form. F25 is the
/// mother claim of F57's prefactor 1.080088 = 2 / 1.851701. Two typed parent
/// edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_{−1} = 4</c>
///         (decay rate), <c>a_0 = 2</c> (dCΨ/dt coefficient), <c>a_3 = 1/4</c>
///         (CΨ crossing threshold).</item>
///   <item><see cref="QuarterAsBilinearMaxvalClaim"/>: typed parent grounding
///         CrossingThreshold = 1/4 on the bilinear-apex maxval axis. Same
///         anchor as F57, Dicke, F60 fold, F62 fold. Added 2026-05-16.</item>
/// </list>
///
/// <para>Tier consistency: F25 is Tier 1 proven (PROOF_MONOTONICITY_CPSI);
/// O(1) evaluation. All three claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>
/// + <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> (which registers
/// QuarterAsBilinearMaxvalClaim) in the builder pipeline.</para></summary>
public static class F25CPsiBellPlusPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF25CPsiBellPlusPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F25CPsiBellPlusPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var quarter = b.Get<QuarterAsBilinearMaxvalClaim>();
            return new F25CPsiBellPlusPi2Inheritance(ladder, quarter);
        });
}
