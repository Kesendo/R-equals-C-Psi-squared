using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class PolarityCubeMapTests
{
    private sealed class FakeBitBClaim : Claim, IZ2AxisClaim
    {
        public FakeBitBClaim() : base("fake BitB", Tier.Tier1Derived, "test") { }
        public override string DisplayName => "fake BitB";
        public override string Summary => "fake BitB";
        public Z2Axis Z2Axis => Z2Axis.BitB;
        public Claim? BitATwin => null;
    }

    private sealed class FakeBitAClaim : Claim, IZ2AxisClaim
    {
        public FakeBitAClaim() : base("fake BitA", Tier.Tier1Derived, "test") { }
        public override string DisplayName => "fake BitA";
        public override string Summary => "fake BitA";
        public Z2Axis Z2Axis => Z2Axis.BitA;
        public Claim? BitATwin => null;
    }

    private sealed class FakeKlein2Claim : Claim, IZ2AxisClaim
    {
        public FakeKlein2Claim() : base("fake Klein2", Tier.Tier1Derived, "test") { }
        public override string DisplayName => "fake Klein2";
        public override string Summary => "fake Klein2";
        public Z2Axis Z2Axis => Z2Axis.Klein2;
        public Claim? BitATwin => null;
    }

    [Fact]
    public void CountsByAxis_AggregatedCorrectly()
    {
        var claims = new List<IZ2AxisClaim>
        {
            new FakeBitBClaim(),
            new FakeBitBClaim(),
            new FakeBitBClaim(),
            new FakeBitAClaim(),
            new FakeKlein2Claim(),
        };

        var map = new PolarityCubeMap(claims);

        Assert.Equal(3, map.BitBClaims.Count);
        Assert.Single(map.BitAClaims);
        Assert.Single(map.Klein2Claims);
        Assert.Empty(map.YParityClaims);
        Assert.Empty(map.Cubic3Claims);
        Assert.Empty(map.NotApplicableClaims);
        Assert.Equal(5, map.TotalClaims);
    }

    [Fact]
    public void OpenBitATwinSlots_CountsBitBClaimsWithNullTwin()
    {
        var claims = new List<IZ2AxisClaim>
        {
            new FakeBitBClaim(),
            new FakeBitBClaim(),
            new FakeBitAClaim(),
        };

        var map = new PolarityCubeMap(claims);

        Assert.Equal(2, map.OpenBitATwinSlots);
        Assert.Equal(0, map.FilledBitATwinSlots);
        Assert.Equal(0.0, map.TwinCoverageRatio);
    }

    [Fact]
    public void EmptyClaimList_GivesEmptyMap()
    {
        var map = new PolarityCubeMap(new List<IZ2AxisClaim>());

        Assert.Equal(0, map.TotalClaims);
        Assert.Equal(0, map.OpenBitATwinSlots);
        Assert.Equal(0.0, map.TwinCoverageRatio);
        Assert.Empty(map.UnfilledTwinSlotNames);
    }

    [Fact]
    public void Constructor_ThrowsOnNull()
    {
        Assert.Throws<ArgumentNullException>(() => new PolarityCubeMap(null!));
    }

    [Fact]
    public void RealPi2KnowledgeBase_ConstructsPolarityCubeMap()
    {
        // Smoke test with the real Pi2KnowledgeBase.
        var chain = new RCPsiSquared.Core.ChainSystems.ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var kb = new Pi2KnowledgeBase(chain);

        Assert.NotNull(kb.PolarityCubeMap);
        Assert.True(kb.PolarityCubeMap.TotalClaims >= 1);
    }
}
