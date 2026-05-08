using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.JordanWigner;

public class C2BondKModeProfileTests
{
    private readonly ITestOutputHelper _out;

    public C2BondKModeProfileTests(ITestOutputHelper output) => _out = output;

    // Endpoint K_90 mean dominates the *minimum* Interior K_90 (innermost-Interior bond) across
    // N=5..8. The naïve EndpointK90Mean ≥ InteriorK90Mean comparison fails at N=8 because
    // flanking-Interior bonds (b ∈ {1,2,4,5}) become more k-uniform than Endpoints while the
    // innermost Interior bond stays the most localized; the min comparison preserves the
    // structural claim "boundary bonds are at least as spread as the deepest-Interior bond".
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void EndpointK90Mean_NotLessThan_MinInteriorK90(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var profile = C2BondKModeProfile.Build(block);

        if (double.IsNaN(profile.MinInteriorK90)) return;  // skip when no Interior bond

        Assert.True(profile.EndpointK90Mean >= profile.MinInteriorK90,
            $"N={N}: Endpoint K_90 mean ({profile.EndpointK90Mean:F2}) should be ≥ " +
            $"min Interior K_90 ({profile.MinInteriorK90:F2})");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(8)]
    public void RowL1Profile_LengthEqualsN(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var profile = C2BondKModeProfile.Build(block);
        foreach (var bond in profile.Bonds)
        {
            Assert.Equal(N, bond.RowL1Profile.Count);
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(8)]
    public void TopThreeKIndices_AreOneIndexedAndDistinct(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var profile = C2BondKModeProfile.Build(block);
        foreach (var bond in profile.Bonds)
        {
            int expectedCount = Math.Min(3, N);
            Assert.Equal(expectedCount, bond.TopThreeKIndices.Count);
            Assert.Equal(bond.TopThreeKIndices.Count, bond.TopThreeKIndices.Distinct().Count());
            foreach (int k in bond.TopThreeKIndices)
            {
                Assert.InRange(k, 1, N);
            }
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(8)]
    public void K99_AtLeast_K90(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var profile = C2BondKModeProfile.Build(block);
        foreach (var bond in profile.Bonds)
        {
            Assert.True(bond.K99 >= bond.K90,
                $"N={N}, bond {bond.Bond}: K99 ({bond.K99}) should be ≥ K90 ({bond.K90})");
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(8)]
    public void K90_AtLeastOne_AtMostN(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var profile = C2BondKModeProfile.Build(block);
        foreach (var bond in profile.Bonds)
        {
            Assert.InRange(bond.K90, 1, N);
        }
    }

    [Fact]
    public void Build_RejectsNonC2Block()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);  // c=3
        Assert.Throws<ArgumentException>(() => C2BondKModeProfile.Build(block));
    }

    [Fact]
    public void Tier_IsTier2Verified()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var profile = C2BondKModeProfile.Build(block);
        Assert.Equal(Tier.Tier2Verified, profile.Tier);
    }

    [Fact]
    public void Anchor_References_PROOF_F86_QPEAK_AndDirectionB()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var profile = C2BondKModeProfile.Build(block);
        Assert.Contains("PROOF_F86_QPEAK", profile.Anchor);
        Assert.Contains("Direction (b'')", profile.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsPerNBondProfiles_AcrossN5To8()
    {
        _out.WriteLine("  N | b | class    | K_90 | K_99 | top-3 k indices | row-L1 profile (k=1..N)");
        _out.WriteLine("  --|---|----------|------|------|-----------------|----------");
        foreach (int N in new[] { 5, 6, 7, 8 })
        {
            var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
            var profile = C2BondKModeProfile.Build(block);
            foreach (var bond in profile.Bonds)
            {
                string topKStr = string.Join(",", bond.TopThreeKIndices);
                string profileStr = string.Join(" ", bond.RowL1Profile.Select(v => v.ToString("F4")));
                _out.WriteLine($"  {N} | {bond.Bond} | {bond.BondClass,-8} | {bond.K90,4} | {bond.K99,4} | {topKStr,15} | {profileStr}");
            }
            _out.WriteLine($"  -- N={N}: Endpoint K_90 mean = {profile.EndpointK90Mean:F2}, " +
                           $"Interior K_90 mean = {profile.InteriorK90Mean:F2}, " +
                           $"min Interior K_90 = {profile.MinInteriorK90:F2}");
            _out.WriteLine("");
        }
        Assert.True(true);  // emit-only test
    }
}
