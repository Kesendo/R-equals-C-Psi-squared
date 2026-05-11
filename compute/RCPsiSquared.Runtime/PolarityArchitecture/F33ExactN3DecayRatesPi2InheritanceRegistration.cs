using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F33ExactN3DecayRatesPi2Inheritance"/>:
/// F33 N=3 exact rational decay rates rate_1=2γ, rate_2=8γ/3, rate_3=10γ/3.
/// Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c>
///         (WeightOneRateCoefficient in rate_1 = 2γ, same anchor as F1
///         TwoFactor, F50 DecayRateFactor).</item>
///   <item><see cref="F50WeightOneDegeneracyPi2Inheritance"/>: F33's rate_1
///         IS F50's universal weight-1 eigenvalue position (= −2γ),
///         specialized to N=3. F33 → F50 typed parent makes the
///         specialization explicit.</item>
/// </list>
///
/// <para>F33 documents the Absorption Theorem α = 2γ·⟨n_XY⟩ at three
/// rational ⟨n_XY⟩ values: 1 (pure w=1), 4/3 (mix), 5/3 (mix). At N ≥ 4 the
/// internal rates lose their rational closed form (only F50's 2γ and
/// F43's 2(N−1)γ boundary rates remain universal); F33's rationality is
/// N=3-specific.</para>
///
/// <para>Tier consistency: F33 is Tier 1 exact-rational (closed form via
/// N=3 diagonalization); both parent claims Tier1Derived.</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F50WeightOneDegeneracyPi2InheritanceRegistration.RegisterF50WeightOneDegeneracyPi2Inheritance"/>.</para></summary>
public static class F33ExactN3DecayRatesPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF33ExactN3DecayRatesPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F33ExactN3DecayRatesPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f50 = b.Get<F50WeightOneDegeneracyPi2Inheritance>();
            _ = b.Get<AbsorptionTheoremClaim>();
            return new F33ExactN3DecayRatesPi2Inheritance(ladder, f50);
        });
}
