using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.SymmetryFamily;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F93DetuningAntiPalindromicSpectralInvariance"/>;
/// three typed parents (<see cref="JointPopcountSectors"/>, <see cref="F71MirrorBlockRefinement"/>,
/// <see cref="SymmetryFamilyInventory"/>).</summary>
public static class F93DetuningAntiPalindromicSpectralInvarianceRegistration
{
    public static ClaimRegistryBuilder RegisterF93DetuningAntiPalindromicSpectralInvariance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F93DetuningAntiPalindromicSpectralInvariance>(b =>
        {
            var sectors = b.Get<JointPopcountSectors>();
            var f71 = b.Get<F71MirrorBlockRefinement>();
            var inventory = b.Get<SymmetryFamilyInventory>();
            return new F93DetuningAntiPalindromicSpectralInvariance(sectors, f71, inventory);
        });
}
