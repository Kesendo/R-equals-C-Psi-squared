using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F55UniversalAbsorptionDosePi2Inheritance"/>:
/// F55's absorption dose K_death = ln(10), valid above the D6 coupling
/// threshold Q*_gap(N). Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c>
///         (RateMinCoefficient in rate_min = 2γ).</item>
///   <item><see cref="F50WeightOneDegeneracyPi2Inheritance"/>: F55's rate_min
///         IS F50's weight-1 eigenvalue position (= −2γ). The "2" in
///         rate_min and the "2" in 2·ln(10) cancel exactly to give
///         K_death = ln(10). Note F50 supplies EXISTENCE of modes at 2γ;
///         that 2γ is the MINIMUM is D6, and only above Q*_gap(N).</item>
/// </list>
///
/// <para>Tier consistency: F55 is Tier 1 proven from D6 + 99%-absorption
/// convention; both parent claims Tier1Derived.</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F50WeightOneDegeneracyPi2InheritanceRegistration.RegisterF50WeightOneDegeneracyPi2Inheritance"/>.</para></summary>
public static class F55UniversalAbsorptionDosePi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF55UniversalAbsorptionDosePi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F55UniversalAbsorptionDosePi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f50 = b.Get<F50WeightOneDegeneracyPi2Inheritance>();
            _ = b.Get<AbsorptionTheoremClaim>();
            return new F55UniversalAbsorptionDosePi2Inheritance(ladder, f50);
        });
}
