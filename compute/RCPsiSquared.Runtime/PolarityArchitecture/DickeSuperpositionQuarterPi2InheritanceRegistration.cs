using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="DickeSuperpositionQuarterPi2Inheritance"/>:
/// Theorem 1 + 2 of PROOF_BLOCK_CPSI_QUARTER as Pi2-Foundation inheritance.
/// Three parent edges:
///
/// <list type="bullet">
///   <item><see cref="QuarterAsBilinearMaxvalClaim"/>: provides the 1/4 ceiling
///         (= <c>a_3</c>).</item>
///   <item><see cref="HalfAsStructuralFixedPointClaim"/>: provides the 1/2 sector
///         balance at AM-GM saturation (= <c>a_2</c>).</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides both ladder terms via
///         <c>Term(2)</c> and <c>Term(3)</c>.</item>
/// </list>
///
/// <para>This is one of the genuinely Tier 1 derived F86 closed forms — Theorems
/// 1 + 2 are algebraic identities proven directly (Cauchy-Schwarz + AM-GM +
/// Dicke amplitude uniformity). Distinct from the F86 closed-form territory where
/// values were brute-forced or guessed (per Tom 2026-05-09 caveat).</para>
///
/// <para>Tier consistency: all four Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> in the
/// builder pipeline. Both QuarterAsBilinearMaxval and HalfAsStructuralFixedPoint
/// are part of Pi2Family.</para></summary>
public static class DickeSuperpositionQuarterPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterDickeSuperpositionQuarterPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<DickeSuperpositionQuarterPi2Inheritance>(b =>
        {
            var quarter = b.Get<QuarterAsBilinearMaxvalClaim>();
            var half = b.Get<HalfAsStructuralFixedPointClaim>();
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            return new DickeSuperpositionQuarterPi2Inheritance(ladder, quarter, half);
        });
}
