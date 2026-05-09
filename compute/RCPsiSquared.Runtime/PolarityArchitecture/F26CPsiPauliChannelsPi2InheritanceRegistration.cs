using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F26CPsiPauliChannelsPi2Inheritance"/>:
/// F26 generalizes F25 to general Pauli channels (γ_x, γ_y, γ_z). Two parent
/// edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_{−1} = 4</c>
///         (decay rate coefficient).</item>
///   <item><see cref="F25CPsiBellPlusPi2Inheritance"/>: F25 IS F26's
///         γ_x = γ_y = 0 special case (Bell+ Z-dephasing); typed
///         mother-corollary edge. Pattern parallel to F75 → F77, F81 → F82,
///         F25 → F57.</item>
/// </list>
///
/// <para>Tier consistency: F26 is Tier 1 proven (PROOF_MONOTONICITY_CPSI);
/// O(1) evaluation. Both claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F25CPsiBellPlusPi2InheritanceRegistration.RegisterF25CPsiBellPlusPi2Inheritance"/>.</para></summary>
public static class F26CPsiPauliChannelsPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF26CPsiPauliChannelsPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F26CPsiPauliChannelsPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f25 = b.Get<F25CPsiBellPlusPi2Inheritance>();
            return new F26CPsiPauliChannelsPi2Inheritance(ladder, f25);
        });
}
