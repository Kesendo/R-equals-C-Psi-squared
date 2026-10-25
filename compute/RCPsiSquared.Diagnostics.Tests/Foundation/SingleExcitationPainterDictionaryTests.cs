using System;
using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.States;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The α-dictionary in the single-excitation block against the painters' 4^N pipeline, and what it
/// reads past N = 6. Two error models: (1) the α fit stops at a golden-section width of 1e-7 in α, so a letter
/// read twice from trajectories that agree to rounding differs by at most 1e-7/δJ = 5e-6 in f; (2) the fit's
/// 512-point grid seed picks a basin by argmin, and a near-tie between seeds can flip the basin and move f by
/// about 1. The handshake gate 1e-4 sits between the two, so it passes the first and catches the second.</summary>
public class SingleExcitationPainterDictionaryTests
{
    private const double BasinFlipGate = 1e-4;
    // two fits of mirror-related inputs each stop within 1e-7 of their minimiser in α, so their letters
    // differ by at most 2·1e-7/δJ = 1e-5 in f where the misfit has a well-defined minimum
    private const double TwoStops = 2e-7 / SingleExcitationPainterDictionary.DefaultDeltaJ;

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void Letters_AgreeWithTheFullPipeline(int n)
    {
        var block = new SingleExcitationPainterDictionary(n);
        for (int b = 0; b < n - 1; b++)
        {
            var s = new Symphony(n: n, j: 1.0, gamma: 0.05, hType: HamiltonianType.XY,
                                 initialState: InitialStateKind.BondingMode, defectBond: b, deltaJ: 0.02);
            var pm = ((IInspectable)s).Children.OfType<PaintersMovement>().Single();
            for (int i = 0; i < n; i++)
            {
                Assert.True(Math.Abs(pm.F[i] - block.F[b][i]) < BasinFlipGate,
                    $"N={n} bond {b} site {i}: pipeline f = {pm.F[i]}, block f = {block.F[b][i]}");
                Assert.Equal(pm.Reliable[i], block.Reliable[b][i]);
            }
        }
    }

    /// <summary>The premise of the chain-mirror identity: the reliability flags are mirror-consistent exactly,
    /// and on the reliable sites the letter of bond N − 2 − b is R applied to the letter of bond b to two
    /// stopping widths (on a rejected site the misfit is flat and rounding moves its minimiser, which is why
    /// the painters reject it). At even N the centre bond is its own mirror, so its letter is mirror-even to the
    /// same bound.</summary>
    [Theory]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    [InlineData(10)]
    [InlineData(11)]
    [InlineData(12)]
    public void MirrorBondsLetter_IsTheMirrorImage(int n)
    {
        var d = new SingleExcitationPainterDictionary(n);
        for (int b = 0; b < n - 1; b++)
            for (int i = 0; i < n; i++)
            {
                Assert.Equal(d.Reliable[b][n - 1 - i], d.Reliable[n - 2 - b][i]);
                if (!d.Reliable[n - 2 - b][i]) continue;
                Assert.True(Math.Abs(d.F[n - 2 - b][i] - d.F[b][n - 1 - i]) <= TwoStops,
                    $"N={n} bond {b} site {i}: {d.F[n - 2 - b][i]} against {d.F[b][n - 1 - i]}");
            }
    }

    /// <summary>The identity itself is algebra: for any f and its mirror image Rf, cos(f, Rf) = 1 − 2w(f).
    /// Checked on each letter and its exact mirror image, to machine rounding.</summary>
    [Theory]
    [InlineData(7)]
    [InlineData(12)]
    public void CosOfALetterWithItsMirrorImage_IsOneMinusTwiceTheOddWeight(int n)
    {
        var d = new SingleExcitationPainterDictionary(n);
        for (int b = 0; b < n - 1; b++)
        {
            var f = d.F[b];
            var rf = Enumerable.Range(0, n).Select(i => f[n - 1 - i]).ToArray();
            double cos = SingleExcitationPainterDictionary.Cos(f, rf);
            Assert.True(Math.Abs(cos - (1 - 2 * d.MirrorOddWeight(b))) < 1e-12,
                $"N={n} bond {b}: cos = {cos}, 1 − 2w = {1 - 2 * d.MirrorOddWeight(b)}");
        }
    }

    /// <summary>The reading over the sites reliable for every bond, pinned to the printed three decimals
    /// (the fit's width moves a cosine by far less than the 5e-4 of rounding), with the class of the worst
    /// pair: a chain-mirror pair (a + b = N − 2) or a pair across the centre (a &lt; (N − 2)/2 &lt; b).</summary>
    [Theory]
    [InlineData(3, -0.976, true)]
    [InlineData(5, -0.965, false)]
    [InlineData(6, -0.962, true)]
    [InlineData(7, -0.940, false)]
    [InlineData(8, -0.927, false)]
    [InlineData(9, -0.992, true)]
    [InlineData(10, -0.950, false)]
    [InlineData(11, -0.914, true)]
    [InlineData(12, -0.961, true)]
    public void OverReliableSites_TheWorstPair(int n, double expected, bool mirrorPair)
    {
        var w = new SingleExcitationPainterDictionary(n).WorstAntiCollinearPair(reliableOnly: true);
        Assert.True(Math.Abs(w.Cos - expected) < 5e-4, $"N={n}: {w.Cos:F4} at ({w.A},{w.B})");
        if (mirrorPair) Assert.Equal(n - 2, w.A + w.B);
        else Assert.True(w.A < (n - 2) / 2.0 && (n - 2) / 2.0 < w.B, $"N={n}: ({w.A},{w.B}) is not across the centre");
    }

    /// <summary>N = 4: every site reliable, the one chain-mirror pair (0,2) at w = 0.707, cos = −0.414, and the
    /// worst pair (0,1) at −0.541, an edge bond with the mirror-even centre (neither class above).</summary>
    [Fact]
    public void AtNFour_TheThreeBondInventory()
    {
        var d = new SingleExcitationPainterDictionary(4);
        Assert.All(d.SitesReliableForEveryBond(), Assert.True);
        Assert.True(Math.Abs(d.MirrorOddWeight(0) - 0.707) < 5e-4, $"w(0) = {d.MirrorOddWeight(0)}");
        Assert.True(Math.Abs(SingleExcitationPainterDictionary.Cos(d.F[0], d.F[2]) - (-0.414)) < 5e-4);
        var w = d.WorstAntiCollinearPair(reliableOnly: true);
        Assert.True(Math.Abs(w.Cos - (-0.541)) < 5e-4 && (w.A, w.B) == (0, 1), $"{w.Cos:F4} at ({w.A},{w.B})");
    }

    /// <summary>Over all sites N = 6 reads −0.378 and N = 10 −0.701; the sites the painters' rule rejects carry
    /// those values (two of six at N = 6, with |f| up to 216; two of ten at N = 10, |f| about 42).</summary>
    [Theory]
    [InlineData(6, -0.378, 4, 215.6)]
    [InlineData(10, -0.701, 8, 41.8)]
    public void OverAllSites_TheWeakReadingIsCarriedByRejectedSites(int n, double allSites, int reliableSites, double maxF)
    {
        var d = new SingleExcitationPainterDictionary(n);
        Assert.True(Math.Abs(d.WorstAntiCollinearPair(reliableOnly: false).Cos - allSites) < 5e-4);
        var keep = d.SitesReliableForEveryBond();
        Assert.Equal(reliableSites, keep.Count(x => x));
        double maxRejected = Enumerable.Range(0, n).Where(i => !keep[i])
            .Max(i => Enumerable.Range(0, n - 1).Max(b => Math.Abs(d.F[b][i])));
        Assert.True(Math.Abs(maxRejected - maxF) < 0.05, $"N={n}: max |f| on rejected sites {maxRejected:F2}");
    }

    /// <summary>N = 5: the painted chain-mirror pairs read −0.783 and −0.894; the worst pair, −0.965, is
    /// (0,2)/(1,3), across the centre.</summary>
    [Fact]
    public void AtNFive_TheMirrorPairsAreNotTheWorstPair()
    {
        var d = new SingleExcitationPainterDictionary(5);
        Assert.True(Math.Abs(SingleExcitationPainterDictionary.Cos(d.F[0], d.F[3]) - (-0.783)) < 5e-4);
        Assert.True(Math.Abs(SingleExcitationPainterDictionary.Cos(d.F[1], d.F[2]) - (-0.894)) < 5e-4);
        Assert.True(Math.Abs(SingleExcitationPainterDictionary.Cos(d.F[1], d.F[3]) - (-0.965)) < 5e-4);
    }

    [Fact]
    public void AtNTwelve_BondThreeIsNearlyMirrorOdd() =>
        Assert.True(Math.Abs(new SingleExcitationPainterDictionary(12).MirrorOddWeight(3) - 0.980) < 5e-4);

    [Fact]
    public void Guard_RejectsOutOfRange()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new SingleExcitationPainterDictionary(2));
        Assert.Throws<ArgumentOutOfRangeException>(() => new SingleExcitationPainterDictionary(25));
        Assert.Throws<ArgumentOutOfRangeException>(() => new SingleExcitationPainterDictionary(5, deltaJ: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new SingleExcitationPainterDictionary(5, timePoints: 1));
    }
}
