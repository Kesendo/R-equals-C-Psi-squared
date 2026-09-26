using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>F154 inherits its per-site disagreement rate from Absorption and its
/// mirror-balanced profile locus from F91. Its (1,1) half stands on the joint-popcount grading,
/// on F143's Gram (the zero-frequency room's spectrum) and on F140's frozen root (the lower
/// endpoint at every J). None of the five supplies the interval theorem or the N=11 obstruction.
/// Resolution is topological, so the registration order among siblings is free.</summary>
public static class CompressedDensityLocusClaimRegistration
{
    public static ClaimRegistryBuilder RegisterCompressedDensityLocusClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<CompressedDensityLocusClaim>(context => new CompressedDensityLocusClaim(
            context.Get<AbsorptionTheoremClaim>(),
            context.Get<F71AntiPalindromicGammaSpectralInvariance>(),
            context.Get<JointPopcountSectors>(),
            context.Get<SeedRungGramClaim>(),
            context.Get<FrozenDivisorClaim>()));
}
