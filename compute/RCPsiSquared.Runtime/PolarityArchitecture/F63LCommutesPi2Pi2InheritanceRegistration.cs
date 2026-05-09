using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F63LCommutesPi2Pi2Inheritance"/>:
/// F63's <c>[L, Π²] = 0</c> conservation as F38 sister-claim. Two parent
/// edges:
///
/// <list type="bullet">
///   <item><see cref="F38Pi2InvolutionPi2Inheritance"/>: F63 says L respects
///         what F38 defines; the per-block dimension agrees with F38's
///         eigenspace dimension divided by 2 (Π² eigenspace splits each bit_a
///         sector). Higher inheritance edge than other F-formulas.</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_{−1} = 4</c>
///         (BlockCount) and <c>a_{3−2N} = 4^(N−1)</c> (PerBlockDimension).</item>
/// </list>
///
/// <para>Tier consistency: F63 is Tier 1 proven analytically (six-line proof
/// in PROOF_BIT_B_PARITY_SYMMETRY); commutator vanishes identically. All three
/// claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F88PopcountCoherenceRegistration.RegisterF88PopcountCoherence"/> +
/// <see cref="F88StaticDyadicAnchorRegistration.RegisterF88StaticDyadicAnchor"/> +
/// <see cref="Pi2OperatorSpaceMirrorRegistration.RegisterPi2OperatorSpaceMirror"/> +
/// <see cref="Pi2I4MemoryLoopRegistration.RegisterPi2I4MemoryLoop"/> +
/// <see cref="F38Pi2InvolutionPi2InheritanceRegistration.RegisterF38Pi2InvolutionPi2Inheritance"/>
/// in the builder pipeline.</para></summary>
public static class F63LCommutesPi2Pi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF63LCommutesPi2Pi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F63LCommutesPi2Pi2Inheritance>(b =>
        {
            var f38 = b.Get<F38Pi2InvolutionPi2Inheritance>();
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            return new F63LCommutesPi2Pi2Inheritance(f38, ladder);
        });
}
