using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F3DecayRateBoundsPi2Inheritance"/>:
/// F3 decay rate bounds (min = 2γ, max = 2(N−1)γ, bandwidth = 2(N−2)γ,
/// XOR = 2Nγ). Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c>
///         (RateCoefficient — universal "2" in all four bounds).</item>
///   <item><see cref="F50WeightOneDegeneracyPi2Inheritance"/>: F3's min
///         rate IS F50's universal weight-1 eigenvalue position (= −2γ).</item>
/// </list>
///
/// <para>F3 is a Tier 1 corollary of the Absorption Theorem α = 2γ·⟨n_XY⟩.
/// All four bounds (min, max, bandwidth, XOR) read 2γ multiplied by the
/// appropriate weight count (1, N−1, N−2, N).</para>
///
/// <para>Tier consistency: F3 Tier 1 derived; both parent claims Tier1Derived.</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F50WeightOneDegeneracyPi2InheritanceRegistration.RegisterF50WeightOneDegeneracyPi2Inheritance"/>.</para></summary>
public static class F3DecayRateBoundsPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF3DecayRateBoundsPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F3DecayRateBoundsPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f50 = b.Get<F50WeightOneDegeneracyPi2Inheritance>();
            return new F3DecayRateBoundsPi2Inheritance(ladder, f50);
        });
}
