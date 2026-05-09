using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F73SpatialSumPurityClosurePi2Inheritance"/>:
/// F73's spatial-sum coherence purity closure. F73 sits in F70 family with F72.
/// Three parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_2 = 1/2</c>
///         (baseline), <c>a_{−1} = 4</c> (decay rate), <c>a_0 = 2</c>
///         (per-site coefficient).</item>
///   <item><see cref="F70DeltaNSelectionRulePi2Inheritance"/>: cited
///         scaffolding. F73's (vac, SE) block focus IS F70's |ΔN| ≤ 1
///         single-site bound applied to coherent probes.</item>
///   <item><see cref="F72BlockDiagonalPurityPi2Inheritance"/>: cited
///         scaffolding. F73's spatial-sum is the CC contribution F72
///         isolates as the (vac, SE) coherence block.</item>
/// </list>
///
/// <para>Tier consistency: F73 is Tier 1 proven (general U(1) case +
/// alternative XY derivation); verified ~10⁻¹⁶ deviation. All four claims
/// Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F70DeltaNSelectionRulePi2InheritanceRegistration.RegisterF70DeltaNSelectionRulePi2Inheritance"/> +
/// <see cref="F72BlockDiagonalPurityPi2InheritanceRegistration.RegisterF72BlockDiagonalPurityPi2Inheritance"/>.</para></summary>
public static class F73SpatialSumPurityClosurePi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF73SpatialSumPurityClosurePi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F73SpatialSumPurityClosurePi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f70 = b.Get<F70DeltaNSelectionRulePi2Inheritance>();
            var f72 = b.Get<F72BlockDiagonalPurityPi2Inheritance>();
            return new F73SpatialSumPurityClosurePi2Inheritance(ladder, f70, f72);
        });
}
