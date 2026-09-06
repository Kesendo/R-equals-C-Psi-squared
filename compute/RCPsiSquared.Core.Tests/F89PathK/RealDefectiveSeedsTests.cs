using System;
using System.Linq;
using RCPsiSquared.Core.F89PathK;
using Xunit;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The typed registry of (1,2)-block count-change loci used by the shell census, with the
/// independent character-certification boundary pinned explicitly,
/// (gate RealSeedCensusTests / SLOW_SEEDCENSUS; table in experiments/F89_PATH_K_DIABOLIC.md).</summary>
public class RealDefectiveSeedsTests
{
    [Fact]
    public void Counts_MatchTheCountChangeCensus()
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
    public void CharacterCertification_StopsBeforeN11()
    {
        Assert.All(RealDefectiveSeeds.All.Where(s => s.N is 5 or 7 or 9), s => Assert.True(s.CharacterCertified));
        Assert.All(RealDefectiveSeeds.ForN(11), s => Assert.False(s.CharacterCertified));
    }

    [Fact]
    public void NewLocus_IsNotCharacterCertifiedUnlessExplicitlyOptedIn()
    {
        var futureCandidate = new RealSeed(13, 1.0, -4.0, +1, "test candidate");
        Assert.False(futureCandidate.CharacterCertified);
    }

    [Fact]
    public void EveryRecordedLocus_HasRealLambdaInTheOneTwoWindow()
    {
        // the (1,2) block's Bendixson window is [-6,-2] for every N (n_diff in {1,3});
        // the window-edge lemma sharpens it to the OPEN interval for certified seeds;
        // every recorded count-change locus also lies inside that same numerical window.
        Assert.All(RealDefectiveSeeds.All, s =>
        {
            Assert.True(s.LambdaA > -6.0 && s.LambdaA < -2.0);
            Assert.True(s.QStar > 0.2 && s.QStar < 3.0);   // the census window
            Assert.True(s.RParity is +1 or -1);
        });
    }
}
