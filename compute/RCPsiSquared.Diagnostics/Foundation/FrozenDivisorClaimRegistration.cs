using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="FrozenDivisorClaim"/> (F140; Tier1Derived). Two typed parents:
/// <see cref="F71AntiPalindromicGammaSpectralInvariance"/>, the anti-palindromic γ-locus the
/// divisor lives on, and the reading that the whole R90 orbit shares one diagonal-block spectrum
/// with the uniform profile; and <see cref="JointPopcountSectors"/>, because "the corner block"
/// presupposes the exact (N+1)² block decomposition, which is also what carries the divisor past
/// the N = 8 spectral wall. Resolution is topological, so registration order among siblings is
/// immaterial.</summary>
public static class FrozenDivisorClaimRegistration
{
    public static ClaimRegistryBuilder RegisterFrozenDivisorClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<FrozenDivisorClaim>(b =>
            new FrozenDivisorClaim(
                b.Get<F71AntiPalindromicGammaSpectralInvariance>(),
                b.Get<JointPopcountSectors>()));
}
