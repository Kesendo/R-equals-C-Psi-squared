using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F96BornSubdominantSlopesPi2Inheritance"/>:
/// F96 closes the per-outcome Born-deviation table for F94's setup with the
/// subdominant slopes Δ_|01⟩ = Δ_|10⟩ = −(16/9)·K and Δ_|11⟩ = −(8/3)·K, both
/// simple algebraic expressions in F94's 4/3 anchor.
///
/// <list type="bullet">
///   <item><see cref="F94BornDeviationFourThirdsPi2Inheritance"/>: the 4/3 unit
///         that F96 elaborates. F96 expresses both subdominant slopes as clean
///         algebraic functions of F94's coefficient: <c>−(4/3)²</c> for the
///         singly-flipped outcomes and <c>−2·(4/3)</c> for the doubly-flipped.</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: the dyadic ladder providing the
///         <c>a_{−1} = 4</c> structural unit that is inherited through F94 into
///         F96 (via <c>U_2_SingleFlipped_TimesFour</c> and similar implicit
///         denominators in F96's bit-exact integer storage).</item>
/// </list>
///
/// <para>Tier consistency: F96 Tier 1 derived (bit-exact symbolic via Dyson
/// matrix elements M_3, M_5 and unitary matrix elements U_2, U_4). Both typed
/// parents are Tier 1 derived. Tier composition well-formed.</para>
///
/// <para>Requires:
/// <see cref="F94BornDeviationFourThirdsPi2InheritanceRegistration.RegisterF94BornDeviationFourThirdsPi2Inheritance"/>
/// + <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>.</para></summary>
public static class F96BornSubdominantSlopesPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF96BornSubdominantSlopesPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F96BornSubdominantSlopesPi2Inheritance>(b =>
        {
            var f94 = b.Get<F94BornDeviationFourThirdsPi2Inheritance>();
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            return new F96BornSubdominantSlopesPi2Inheritance(f94, ladder);
        });
}
