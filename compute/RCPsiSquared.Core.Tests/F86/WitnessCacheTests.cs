using System;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86;

public class WitnessCacheTests
{
    [Fact]
    public void GetOrComputeC2HwhmRatio_ReturnsSameInstance_ForRepeatCalls()
    {
        var cache = new WitnessCache();
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);

        var first = cache.GetOrComputeC2HwhmRatio(block);
        var second = cache.GetOrComputeC2HwhmRatio(block);

        Assert.Same(first, second);
    }

    [Fact]
    public void GetOrComputeC2HwhmRatio_SharesAcrossClaims_AvoidsDoubleQScan()
    {
        // C2UniversalShapeDerivation and PerF71OrbitKTable both consume C2HwhmRatio.
        // When given the same WitnessCache, both must reuse the cached instance so an
        // F86KnowledgeBase inspection pulling both pays the Q-scan once. We verify by
        // reference-equal witness collections: universalShape.Witnesses reads from the
        // C2HwhmRatio.Witnesses property; the cache's C2HwhmRatio is the same object both
        // claims received.
        var cache = new WitnessCache();
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);

        var universalShape = C2UniversalShapeDerivation.Build(block, cache);
        _ = PerF71OrbitKTable.Build(block, cache);
        var cachedHwhmRatio = cache.GetOrComputeC2HwhmRatio(block);

        Assert.Same(cachedHwhmRatio, universalShape.HwhmRatio);
        Assert.Same(cachedHwhmRatio.Witnesses, universalShape.Witnesses);
    }

    [Fact]
    public void GetOrComputeC2HwhmRatio_ThrowsIfNotC2()
    {
        var cache = new WitnessCache();
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);  // c=3

        Assert.Throws<ArgumentException>(() => cache.GetOrComputeC2HwhmRatio(block));
    }
}
