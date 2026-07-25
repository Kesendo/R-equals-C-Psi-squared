using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Schicht-1 wiring tests for <see cref="FrozenDivisorClaim"/> (F140, the R90 frozen
/// divisor), with typed parent edges to
/// <see cref="F71AntiPalindromicGammaSpectralInvariance"/> (the anti-palindromic γ-locus the
/// divisor lives on) and <see cref="JointPopcountSectors"/> (the block decomposition "the corner
/// block" presupposes).</summary>
public class FrozenDivisorClaimRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterJointPopcountSectors()
            .RegisterF71MirrorBlockRefinement()
            .RegisterF71AntiPalindromicGammaSpectralInvariance();

    [Fact]
    public void RegisterFrozenDivisor_AddsClaim_Tier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterFrozenDivisorClaim()
            .Build();
        Assert.True(registry.Contains<FrozenDivisorClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<FrozenDivisorClaim>().Tier);
    }

    [Fact]
    public void RegisterFrozenDivisor_AncestorsContainBothParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterFrozenDivisorClaim()
            .Build();
        var ancestors = registry.AncestorsOf<FrozenDivisorClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F71AntiPalindromicGammaSpectralInvariance), ancestors);
        Assert.Contains(typeof(JointPopcountSectors), ancestors);
    }

    /// <summary>The closed forms the claim carries, checked against the numbers the proof states:
    /// one frozen mode per balanced pair, twice that many mirror-fixed cells, the pair distances
    /// N + 1 − 2c, and the total valuation 2⌊N²/4⌋.</summary>
    [Theory]
    [InlineData(3, 1, 4)]
    [InlineData(4, 2, 8)]
    [InlineData(5, 2, 12)]
    [InlineData(6, 3, 18)]
    [InlineData(7, 3, 24)]
    [InlineData(10, 5, 50)]
    public void ClosedForms_MatchTheProof(int n, int frozen, int totalValuation)
    {
        Assert.Equal(frozen, FrozenDivisorClaim.FrozenMultiplicity(n));
        Assert.Equal(2 * frozen, FrozenDivisorClaim.MirrorFixedCells(n));
        Assert.Equal(totalValuation, FrozenDivisorClaim.TotalValuation(n));

        int sum = 0;
        for (int c = 1; c <= frozen; c++)
        {
            Assert.Equal(n + 1 - 2 * c, FrozenDivisorClaim.PairDistance(n, c));
            sum += FrozenDivisorClaim.PairDistance(n, c);
        }
        Assert.Equal(totalValuation, 2 * sum);
    }

    /// <summary>The two roots, and the one N where the gamma fold sends the root to itself.
    /// </summary>
    [Fact]
    public void TheTwoRoots_CoincideExactlyAtNEqualsFour()
    {
        const double gBar = 0.09;
        Assert.Equal(-4.0 * gBar, FrozenDivisorClaim.FrozenRoot(gBar), 12);
        for (int n = 3; n <= 10; n++)
        {
            bool same = Math.Abs(FrozenDivisorClaim.FoldedRoot(n, gBar) - FrozenDivisorClaim.FrozenRoot(gBar)) < 1e-12;
            Assert.Equal(n == 4, same);
        }
    }
}
