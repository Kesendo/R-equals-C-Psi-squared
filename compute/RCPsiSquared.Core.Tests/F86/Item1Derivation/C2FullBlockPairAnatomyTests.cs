using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class C2FullBlockPairAnatomyTests
{
    private readonly ITestOutputHelper _out;

    public C2FullBlockPairAnatomyTests(ITestOutputHelper output) => _out = output;

    private static CoherenceBlock C2Block(int N) => new(N: N, n: 1, gammaZero: 0.05);

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Build_AcrossN5To8_LandsTier2Verified(int N)
    {
        var anatomy = C2FullBlockPairAnatomy.Build(C2Block(N));
        Assert.Equal(Tier.Tier2Verified, anatomy.Tier);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void K90Pair_IsFinite_AndPositive(int N)
    {
        var anatomy = C2FullBlockPairAnatomy.Build(C2Block(N));
        Assert.True(anatomy.K90Pair >= 1, $"N={N}: K_90_pair must be ≥ 1");
        int totalPairs = anatomy.Block.Basis.MTotal * anatomy.Block.Basis.MTotal;
        Assert.True(anatomy.K90Pair <= totalPairs, $"N={N}: K_90_pair ≤ dim²");
        Assert.True(anatomy.K99Pair >= anatomy.K90Pair, $"N={N}: K_99 ≥ K_90");
    }

    [Fact]
    public void DirectionBPrimePrime_PairLevelReconnaissance_AtQpeak_AcrossN5To8_EmitTable()
    {
        // The decisive question after the diagonal anatomy showed K_90 grows linearly
        // with N: when we look at the FULL Duhamel pair decomposition at Q=Q_peak (where
        // HWHM is actually read), does the K-mass concentrate in a small N-stable set of
        // pairs? Or does the pair view also show linear-in-N growth?
        _out.WriteLine("Direction (b'') pair-level reconnaissance: full Duhamel pair anatomy");
        _out.WriteLine("at (Q = Q_peak ≈ 2.197·Q_EP, t = t_peak = 1/(4γ₀))");
        _out.WriteLine("");
        _out.WriteLine("  N | dim | dim²   | K_90_pair | K_99_pair | top pair (i,j) | top |contrib| | top frac");
        _out.WriteLine("  --|-----|--------|-----------|-----------|----------------|---------------|---------");
        foreach (int N in new[] { 5, 6, 7, 8 })
        {
            var anatomy = C2FullBlockPairAnatomy.Build(C2Block(N));
            int dim = anatomy.Block.Basis.MTotal;
            int totalPairs = dim * dim;
            var top = anatomy.TopPairs[0];
            _out.WriteLine(
                $"  {N} | {dim,3} | {totalPairs,6} | {anatomy.K90Pair,9} | {anatomy.K99Pair,9} | " +
                $"({top.IndexI,3},{top.IndexJ,3})       | {top.AbsContribution,11:G4}   | {top.FractionOfTotal,7:F4}");
        }
        _out.WriteLine("");
        _out.WriteLine("Reading: K_90_pair small + N-stable ⟹ HWHM is structurally a low-rank pair");
        _out.WriteLine("truncation, Tier1Candidate path with explicit truncation error term realistic.");
        _out.WriteLine("K_90_pair growth ⟹ the Duhamel coherent superposition spreads across many");
        _out.WriteLine("pair channels; a clean closed form needs an alternative structural reduction.");
    }

    [Fact]
    public void Build_RejectsNonC2Block()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);
        Assert.Throws<ArgumentException>(() => C2FullBlockPairAnatomy.Build(block));
    }

    [Theory]
    [InlineData(5)]
    public void TopPairs_LimitedToTopPairsToShow(int N)
    {
        var anatomy = C2FullBlockPairAnatomy.Build(C2Block(N));
        Assert.True(anatomy.TopPairs.Count <= C2FullBlockPairAnatomy.TopPairsToShow);
    }

    [Theory]
    [InlineData(5)]
    public void TopPairs_AreSortedDescendingByContribution(int N)
    {
        var anatomy = C2FullBlockPairAnatomy.Build(C2Block(N));
        for (int k = 1; k < anatomy.TopPairs.Count; k++)
            Assert.True(anatomy.TopPairs[k - 1].AbsContribution >= anatomy.TopPairs[k].AbsContribution,
                $"top-pair sort order violated at k={k}");
    }

    [Theory]
    [InlineData(5)]
    public void Children_IncludeK90PairAndTopPairsGroup(int N)
    {
        IInspectable claim = C2FullBlockPairAnatomy.Build(C2Block(N));
        var labels = claim.Children.Select(c => c.DisplayName).ToList();
        Assert.Contains("K_90_pair", labels);
        Assert.Contains(labels, l => l.StartsWith("top ") && l.Contains("pairs"));
    }
}
