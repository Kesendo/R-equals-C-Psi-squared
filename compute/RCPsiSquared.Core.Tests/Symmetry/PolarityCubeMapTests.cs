using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

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
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var kb = new Pi2KnowledgeBase(chain);

        Assert.NotNull(kb.PolarityCubeMap);
        Assert.True(kb.PolarityCubeMap.TotalClaims >= 1);
    }

    private sealed class FakeBitBClaimTrivialNotYetTyped : Claim, IZ2AxisClaim
    {
        public FakeBitBClaimTrivialNotYetTyped() : base("fake BitB trivial", Tier.Tier1Derived, "test") { }
        public override string DisplayName => "fake BitB trivial";
        public override string Summary => "fake BitB trivial";
        public Z2Axis Z2Axis => Z2Axis.BitB;
        public Claim? BitATwin => null;
        public BitATwinClassification BitATwinStatus => BitATwinClassification.TrivialNotYetTyped;
    }

    private sealed class FakeBitBClaimBitBSpecific : Claim, IZ2AxisClaim
    {
        public FakeBitBClaimBitBSpecific() : base("fake BitB specific", Tier.Tier1Derived, "test") { }
        public override string DisplayName => "fake BitB specific";
        public override string Summary => "fake BitB specific";
        public Z2Axis Z2Axis => Z2Axis.BitB;
        public Claim? BitATwin => null;
        public BitATwinClassification BitATwinStatus => BitATwinClassification.BitBSpecific;
    }

    [Fact]
    public void BitATwinStatus_DefaultBehavior_MapsBitBNullToNeedsDerivation()
    {
        // Default IZ2AxisClaim.BitATwinStatus impl: BitB + null BitATwin → NeedsDerivation.
        IZ2AxisClaim claim = new FakeBitBClaim();
        Assert.Equal(BitATwinClassification.NeedsDerivation, claim.BitATwinStatus);
    }

    [Fact]
    public void BitATwinStatus_DefaultBehavior_MapsNonBitBToNotApplicable()
    {
        // Default IZ2AxisClaim.BitATwinStatus impl: non-BitB axis → NotApplicableForThisAxis.
        IZ2AxisClaim bitA = new FakeBitAClaim();
        IZ2AxisClaim klein2 = new FakeKlein2Claim();
        Assert.Equal(BitATwinClassification.NotApplicableForThisAxis, bitA.BitATwinStatus);
        Assert.Equal(BitATwinClassification.NotApplicableForThisAxis, klein2.BitATwinStatus);
    }

    [Fact]
    public void BitATwinStatusBreakdown_AggregatesAcrossBitBClaims()
    {
        // Three BitB Claims with distinct status overrides + one BitA + one Klein2.
        // BitA/Klein2 should NOT count toward any of the breakdown buckets
        // (they're NotApplicableForThisAxis).
        var claims = new List<IZ2AxisClaim>
        {
            new FakeBitBClaim(),                       // default → NeedsDerivation
            new FakeBitBClaimTrivialNotYetTyped(),     // override → TrivialNotYetTyped
            new FakeBitBClaimBitBSpecific(),           // override → BitBSpecific
            new FakeBitAClaim(),                       // → NotApplicableForThisAxis
            new FakeKlein2Claim(),                     // → NotApplicableForThisAxis
        };

        var map = new PolarityCubeMap(claims);

        Assert.Equal(1, map.TrivialNotYetTypedTwinSlots);
        Assert.Equal(1, map.NeedsDerivationTwinSlots);
        Assert.Equal(1, map.BitBSpecificTwinSlots);

        Assert.Single(map.TrivialNotYetTypedTwinSlotNames);
        Assert.Equal("FakeBitBClaimTrivialNotYetTyped", map.TrivialNotYetTypedTwinSlotNames[0]);

        // OpenBitATwinSlots still counts all 3 BitB (none have BitATwin); the
        // breakdown REFINES the open count, doesn't replace it.
        Assert.Equal(3, map.OpenBitATwinSlots);
    }
}
