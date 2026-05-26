using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="LindbladBitBPiBreakMagnitude"/>
/// (F113, Tier1Derived at N=2, 3, 4; Tier1Candidate general N). F113 is the
/// closed-form magnitude for the F112 polarity-asymmetry counterexample, sitting
/// on the same bit_b axis as its parent.
///
/// <para>Typed ctor parent: <see cref="LindbladBitBPiBalance"/> (F112). F112 must
/// be registered BEFORE F113 so the <c>b.Get&lt;...&gt;()</c> resolves; the
/// registration ordering in <c>KnowledgeRegistryFactory.BuildDefault</c> places
/// this call directly after F112 to keep the bit_b-axis Lindblad-polarity cluster
/// contiguous in the factory.</para></summary>
public static class LindbladBitBPiBreakMagnitudeRegistration
{
    public static ClaimRegistryBuilder RegisterLindbladBitBPiBreakMagnitude(
        this ClaimRegistryBuilder builder) =>
        builder.Register<LindbladBitBPiBreakMagnitude>(b =>
            new LindbladBitBPiBreakMagnitude(
                b.Get<LindbladBitBPiBalance>()));
}
