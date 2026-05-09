using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F72BlockDiagonalPurityPi2Inheritance"/>:
/// F72's <c>Tr(ρ_i²) = 1/2 + P_DD + P_CC</c> as Tier 1 corollary of F70.
/// Two parent edges plus one registration discard:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_2 = 1/2</c>
///         (baseline) and <c>a_0 = 2</c> (block count).</item>
///   <item><see cref="F70DeltaNSelectionRulePi2Inheritance"/>: cited as
///         direct foundation. F72's DD ⊕ CC decomposition IS F70's
///         |ΔN| ≤ 1 single-site bound applied to per-site purity.</item>
///   <item><see cref="HalfAsStructuralFixedPointClaim"/>: registration
///         discard. The 1/2 baseline IS the structural three-faces fixed
///         point; transitively reachable but the discard makes the
///         structural reading explicit.</item>
/// </list>
///
/// <para>Tier consistency: F72 is Tier 1 corollary of F70; verified at N=5
/// w-scan with block-diagonal coupling at machine precision. All four
/// claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (registers HalfAsStructuralFixedPointClaim) +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F70DeltaNSelectionRulePi2InheritanceRegistration.RegisterF70DeltaNSelectionRulePi2Inheritance"/>.</para></summary>
public static class F72BlockDiagonalPurityPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF72BlockDiagonalPurityPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F72BlockDiagonalPurityPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f70 = b.Get<F70DeltaNSelectionRulePi2Inheritance>();
            _ = b.Get<HalfAsStructuralFixedPointClaim>();   // 1/2 baseline = structural fixed point
            return new F72BlockDiagonalPurityPi2Inheritance(ladder, f70);
        });
}
