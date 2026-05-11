using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F66PoleModesPi2Inheritance"/>:
/// F66's dissipation interval [0, 2γ₀] for the uniform XY chain with B at
/// endpoint. Two parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c>, the
///         upper-pole γ₀-coefficient = polynomial root d.</item>
///   <item><see cref="QubitDimensionalAnchorClaim"/>: anchors the "2" as the
///         qubit dimensionality d in d² − 2d = 0; cross-pinning of the
///         Pi2DyadicLadder's a_0 entry.</item>
/// </list>
///
/// <para>Tier consistency: F66 is Tier 1 verified analytically + numerically
/// (per ANALYTICAL_FORMULAS, PROOF_ABSORPTION_THEOREM). All three claims
/// Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (registers QubitDimensionalAnchorClaim) +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> in the
/// builder pipeline.</para></summary>
public static class F66PoleModesPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF66PoleModesPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F66PoleModesPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var qubitAnchor = b.Get<QubitDimensionalAnchorClaim>();
            _ = b.Get<AbsorptionTheoremClaim>();
            return new F66PoleModesPi2Inheritance(ladder, qubitAnchor);
        });
}
