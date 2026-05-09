using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F23XorDrainVanishingFractionPi2Inheritance"/>:
/// F23's combinatorial XOR-drain fraction (N+1)/4^N. Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_{−1} = 4</c>
///         (BaseFactor; same anchor as F25/F65/F73/F76 decay rate
///         numerators).</item>
///   <item><see cref="Pi2OperatorSpaceMirrorClaim"/>: provides the typed
///         4^N ↔ 1/4^N operator-space mirror (pinned table N=1..6); F23's
///         denominator IS Pi2OperatorSpaceMirror's OperatorSpace.</item>
/// </list>
///
/// <para>Tier consistency: F23 is Tier 1 combinatorial (proof in
/// <c>experiments/N_INFINITY_PALINDROME.md</c>); valid for any N,
/// Z-dephasing. Both parent claims Tier1Derived.</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> (which registers
/// Pi2OperatorSpaceMirrorClaim).</para></summary>
public static class F23XorDrainVanishingFractionPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF23XorDrainVanishingFractionPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F23XorDrainVanishingFractionPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var mirror = b.Get<Pi2OperatorSpaceMirrorClaim>();
            return new F23XorDrainVanishingFractionPi2Inheritance(ladder, mirror);
        });
}
