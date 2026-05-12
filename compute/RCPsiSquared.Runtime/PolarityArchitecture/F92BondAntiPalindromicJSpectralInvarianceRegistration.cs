using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.SymmetryFamily;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F92BondAntiPalindromicJSpectralInvariance"/>;
/// three typed parents (<see cref="JointPopcountSectors"/>, <see cref="F71MirrorBlockRefinement"/>,
/// <see cref="SymmetryFamilyInventory"/>).</summary>
public static class F92BondAntiPalindromicJSpectralInvarianceRegistration
{
    public static ClaimRegistryBuilder RegisterF92BondAntiPalindromicJSpectralInvariance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F92BondAntiPalindromicJSpectralInvariance>(b =>
        {
            var sectors = b.Get<JointPopcountSectors>();
            var f71 = b.Get<F71MirrorBlockRefinement>();
            var inventory = b.Get<SymmetryFamilyInventory>();
            return new F92BondAntiPalindromicJSpectralInvariance(sectors, f71, inventory);
        });
}
