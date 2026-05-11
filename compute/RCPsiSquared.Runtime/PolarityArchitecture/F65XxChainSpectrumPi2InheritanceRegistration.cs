using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F65XxChainSpectrumPi2Inheritance"/>:
/// F65 single-excitation spectrum source for F75. Two parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_{−1} = 4</c>
///         (numerator coefficient, F25/F73/F76 sibling).</item>
///   <item><see cref="F66PoleModesPi2Inheritance"/>: F65 + F66 sibling pair
///         at the [0, 2γ₀] interval. F66 says the pole positions; F65 says
///         the single-excitation rates fill the interval but never reach
///         the upper pole 2γ₀.</item>
/// </list>
///
/// <para>Tier consistency: F65 is Tier 1 proven; verified N=3..30 to
/// 1.2·10⁻¹⁵. All three claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F66PoleModesPi2InheritanceRegistration.RegisterF66PoleModesPi2Inheritance"/>.</para></summary>
public static class F65XxChainSpectrumPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF65XxChainSpectrumPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F65XxChainSpectrumPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f66 = b.Get<F66PoleModesPi2Inheritance>();
            _ = b.Get<AbsorptionTheoremClaim>();
            return new F65XxChainSpectrumPi2Inheritance(ladder, f66);
        });
}
