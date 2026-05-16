using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F62WStateBornBelowFoldPi2Inheritance"/>:
/// F62's CΨ(0)_W_N closed form. F62 is the W-state companion to F60 (GHZ).
/// Three typed parent edges, all explicit as of 2026-05-16:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_3 = 1/4</c> (fold)
///         and <c>a_0 = 2</c> (numerator coefficient) via live <c>Term(...)</c>
///         reads.</item>
///   <item><see cref="QuarterAsBilinearMaxvalClaim"/>: typed parent grounding
///         the "fold IS the bilinear-apex maxval = (1/2)²" reading; same
///         anchor as F57, Dicke, F60.</item>
///   <item><see cref="F61BitAParityPi2Inheritance"/>: typed parent grounding
///         "Combined with the Parity Selection Rule (F61), this proves that
///         single-excitation states never cross CΨ = 1/4." F61 + F62
///         together close the SE-regime structural-outside-framework reading.</item>
/// </list>
///
/// <para>Tier consistency: F62 is Tier 1 analytical; verified numerically
/// N=2..10. All four claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (registers QuarterAsBilinearMaxvalClaim) +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// the F38/F63/F61 chain (since F61 depends on F63 depends on F38).</para></summary>
public static class F62WStateBornBelowFoldPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF62WStateBornBelowFoldPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F62WStateBornBelowFoldPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var quarter = b.Get<QuarterAsBilinearMaxvalClaim>();
            var f61 = b.Get<F61BitAParityPi2Inheritance>();
            return new F62WStateBornBelowFoldPi2Inheritance(ladder, quarter, f61);
        });
}
