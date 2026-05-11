using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F64CavityModeExposurePi2Inheritance"/>:
/// F64's effective γ from single-site exposure γ_eff = γ_B·|a_B|². Two typed
/// parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c>
///         (EigenvalueConventionCoefficient in α = 2γ_B·|a_B|²; same anchor
///         as F1 TwoFactor, F50 DecayRateFactor, F44 SumCoefficient).</item>
///   <item><see cref="QuarterAsBilinearMaxvalClaim"/>: provides the
///         g(r=1) = 1/4 special value at uniform J for N=3 chain. F64's
///         uniform-J ceiling sits exactly at the bilinear-apex maxval.</item>
/// </list>
///
/// <para>Tier consistency: F64 is Tier 1-2 (analytical closed form +
/// numerically verified at N=3 chain to 1.8%, N=4 chain to 0.0003, and
/// across 5 topologies + non-uniform J at N=5, 7); both parent claims
/// Tier1Derived.</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> (registers
/// QuarterAsBilinearMaxvalClaim).</para></summary>
public static class F64CavityModeExposurePi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF64CavityModeExposurePi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F64CavityModeExposurePi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var quarter = b.Get<QuarterAsBilinearMaxvalClaim>();
            _ = b.Get<AbsorptionTheoremClaim>();
            return new F64CavityModeExposurePi2Inheritance(ladder, quarter);
        });
}
