using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="BlockSpectrumPerformanceWitness"/>: empirical
/// timing witness that the per-block eig path runs at least 2× faster than full-L eig at
/// N=5,6 with multiset-equal spectrum. Has exactly one typed parent —
/// <see cref="LiouvillianBlockSpectrum"/> — whose structural correctness this witness
/// extends with an empirical performance bound.</summary>
public static class BlockSpectrumPerformanceWitnessRegistration
{
    public static ClaimRegistryBuilder RegisterBlockSpectrumPerformanceWitness(
        this ClaimRegistryBuilder builder) =>
        builder.Register<BlockSpectrumPerformanceWitness>(b =>
        {
            var blockSpectrum = b.Get<LiouvillianBlockSpectrum>();
            return new BlockSpectrumPerformanceWitness(blockSpectrum);
        });
}
