using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>F154 inherits its per-site disagreement rate from Absorption and its
/// mirror-balanced profile locus from F91. Neither parent supplies the mixed-space claim.</summary>
public static class CompressedDensityLocusClaimRegistration
{
    public static ClaimRegistryBuilder RegisterCompressedDensityLocusClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<CompressedDensityLocusClaim>(context => new CompressedDensityLocusClaim(
            context.Get<AbsorptionTheoremClaim>(),
            context.Get<F71AntiPalindromicGammaSpectralInvariance>()));
}
