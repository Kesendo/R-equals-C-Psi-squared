using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F86QEpPi2Inheritance"/>: F86 Q_EP's "2"
/// numerator as Pi2-Foundation <c>a_0 = d</c>. Two parent edges:
///
/// <list type="bullet">
///   <item><see cref="QEpLaw"/>: the F86a closed form (parameterised by g_eff).</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c>.</item>
/// </list>
///
/// <para>Tier consistency: all three Tier1Derived.</para>
///
/// <para>Requires: <see cref="F86Main.F86MainRegistration.RegisterF86Main"/> for
/// QEpLaw + <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>.</para></summary>
public static class F86QEpPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF86QEpPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F86QEpPi2Inheritance>(b =>
        {
            var qEp = b.Get<QEpLaw>();
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            return new F86QEpPi2Inheritance(ladder, qEp);
        });
}
