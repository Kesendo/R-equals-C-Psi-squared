using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F61BitAParityPi2Inheritance"/>:
/// F61's <c>[L, Π²_X] = 0</c> bit_a parity selection rule, sister-claim to
/// F63's bit_b parity. Two parent edges:
///
/// <list type="bullet">
///   <item><see cref="F63LCommutesPi2Pi2Inheritance"/>: F61 inherits the
///         shared 4-block structure (BlockCount, PerBlockDimension) from
///         F63's typed claim. The two together complete the C₂ × C₂ pair
///         as standalone Pi2-Inheritance Claims.</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_{−1} = 4</c>
///         and <c>a_{3−2N} = 4^(N−1)</c> directly, paralleling F63.</item>
/// </list>
///
/// <para>Tier consistency: F61 is Tier 1 proven (PROOF_PARITY_SELECTION_RULE);
/// verified 64 configurations N=2..6. All three claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: same pipeline as F63 plus
/// <see cref="F63LCommutesPi2Pi2InheritanceRegistration.RegisterF63LCommutesPi2Pi2Inheritance"/>.</para></summary>
public static class F61BitAParityPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF61BitAParityPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F61BitAParityPi2Inheritance>(b =>
        {
            var f63 = b.Get<F63LCommutesPi2Pi2Inheritance>();
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            return new F61BitAParityPi2Inheritance(f63, ladder);
        });
}
