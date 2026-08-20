using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="PhysicalGeneratorPolarityBreak"/>
/// (F155, Tier1Derived, universal N): the polarity break of the physical (no-jump) generator
/// ρ ↦ −i(Hρ − ρH†) as one diagonal bilinear form in the Pauli basis.
///
/// <para>Typed ctor parent: <see cref="LindbladBitBPiBreakMagnitude"/> (F113), whose closed form
/// is this law's Z-diagonal special case, so F113 must be registered before this call. It sits
/// directly after F113 to keep the bit_b-axis polarity cluster contiguous.</para>
///
/// <para>On the same bit_b axis as its parent, but with a different twin verdict:
/// <see cref="BitATwinClassification.CoveredByHadamardDuality"/> rather than F113's
/// <see cref="BitATwinClassification.BitBSpecific"/>. The two are about different objects, and
/// the reason is on the claim's <c>BitATwinStatus</c>. PolarityCubeMap's covered-by-Hadamard
/// count grows by 1.</para></summary>
public static class PhysicalGeneratorPolarityBreakRegistration
{
    public static ClaimRegistryBuilder RegisterPhysicalGeneratorPolarityBreak(
        this ClaimRegistryBuilder builder) =>
        builder.Register<PhysicalGeneratorPolarityBreak>(b =>
            new PhysicalGeneratorPolarityBreak(
                b.Get<LindbladBitBPiBreakMagnitude>()));
}
