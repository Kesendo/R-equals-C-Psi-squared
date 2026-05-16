using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F86TPeakPi2Inheritance"/>: F86 t_peak's
/// "4" denominator as Pi2-Foundation <c>a_{−1} = d²</c> for N=1. Two parent edges:
///
/// <list type="bullet">
///   <item><see cref="TPeakLaw"/>: the F86a closed form (parameterised by γ₀).</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_{−1} = 4</c>.</item>
/// </list>
///
/// <para>Tier consistency: all three Tier1Derived.</para>
///
/// <para>Requires: <see cref="F86Main.F86MainRegistration.RegisterF86Main"/> for
/// TPeakLaw + <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>.</para></summary>
public static class F86TPeakPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF86TPeakPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F86TPeakPi2Inheritance>(b =>
        {
            var tPeak = b.Get<TPeakLaw>();
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var quarter = b.Get<QuarterAsBilinearMaxvalClaim>();
            return new F86TPeakPi2Inheritance(ladder, tPeak, quarter);
        });
}
