using System;
using System.Linq;
using RCPsiSquared.Core.F89PathK;
using Xunit;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The typed registry of real defective (1,2)-block seeds (the containment corollary's per-N
/// input and the step-3 shell census's probe loci), pinned against the PT-break count-change census
/// (gate RealSeedCensusTests / SLOW_SEEDCENSUS; table in experiments/F89_PATH_K_DIABOLIC.md).</summary>
public class RealDefectiveSeedsTests
{
    [Fact]
    public void Counts_MatchTheSeedCensus()
    {
        Assert.Equal(4, RealDefectiveSeeds.ForN(5).Count());
        Assert.Equal(6, RealDefectiveSeeds.ForN(7).Count());
        Assert.Equal(7, RealDefectiveSeeds.ForN(9).Count());
        Assert.Equal(9, RealDefectiveSeeds.ForN(11).Count());
    }

    [Fact]
    public void N5Anchors_MatchTheGateValues()
    {
        var n5 = RealDefectiveSeeds.ForN(5).ToList();
        Assert.Contains(n5, s => Math.Abs(s.QStar - 0.620878) < 1e-6 && s.RParity == +1);
        Assert.Contains(n5, s => Math.Abs(s.QStar - 1.077615) < 1e-6 && s.RParity == +1);
        Assert.Contains(n5, s => Math.Abs(s.QStar - 2.804888) < 1e-6 && s.RParity == -1);
        Assert.Contains(n5, s => Math.Abs(s.QStar - 0.643037) < 1e-6 && s.RParity == -1);
    }

    [Fact]
    public void EverySeed_HasRealLambdaInTheOneTwoWindow()
    {
        // the (1,2) block's Bendixson window is [-6,-2] for every N (n_diff in {1,3});
        // the window-edge lemma sharpens it to the OPEN interval — every recorded seed obeys it
        Assert.All(RealDefectiveSeeds.All, s =>
        {
            Assert.True(s.LambdaA > -6.0 && s.LambdaA < -2.0);
            Assert.True(s.QStar > 0.2 && s.QStar < 3.0);   // the census window
            Assert.True(s.RParity is +1 or -1);
        });
    }
}
