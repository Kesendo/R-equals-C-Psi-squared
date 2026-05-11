using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F71AntiPalindromicGammaSpectralInvariance"/>: the
/// anti-palindromic-γ spectral-invariance witness. Has TWO typed parents —
/// <see cref="JointPopcountSectors"/> (γ-blind: joint-popcount block-diagonality independent
/// of γ_l) and <see cref="F71MirrorBlockRefinement"/> (γ-symmetric: exact only when
/// γ_l = γ_{N-1-l}, of which anti-palindromy γ_l + γ_{N-1-l} = 2·γ_avg is a strict
/// weakening). Together they bracket the empirical finding that anti-palindromic γ leaves
/// the spectrum invariant despite breaking F71 in operator form.</summary>
public static class F71AntiPalindromicGammaSpectralInvarianceRegistration
{
    public static ClaimRegistryBuilder RegisterF71AntiPalindromicGammaSpectralInvariance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F71AntiPalindromicGammaSpectralInvariance>(b =>
        {
            var sectors = b.Get<JointPopcountSectors>();
            var f71 = b.Get<F71MirrorBlockRefinement>();
            return new F71AntiPalindromicGammaSpectralInvariance(sectors, f71);
        });
}
