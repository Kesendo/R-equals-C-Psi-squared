using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F71MirrorSymmetryPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterHalfIntegerMirror(N: 5);

    [Fact]
    public void RegisterF71_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF71MirrorSymmetryPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F71MirrorSymmetryPi2Inheritance>());
    }

    [Fact]
    public void RegisterF71_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF71MirrorSymmetryPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F71MirrorSymmetryPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF71_AncestorsContainHalfIntegerMirror_FirstDirectEdge()
    {
        // Before F71, HalfIntegerMirrorClaim had 0 descendants. F71 is the first
        // F-formula on the half-integer-mirror axis (per Tom 2026-05-09
        // mirror-map check).
        var registry = BuildBaseRegistry()
            .RegisterF71MirrorSymmetryPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F71MirrorSymmetryPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(HalfIntegerMirrorClaim), ancestors);
    }

    [Fact]
    public void RegisterF71_HalfIntegerMirrorNowHasDescendants()
    {
        // The mirror-map detector showed HalfIntegerMirrorClaim with 0
        // descendants; F71 fills the gap. Cross-registry verification.
        var registry = BuildBaseRegistry()
            .RegisterF71MirrorSymmetryPi2Inheritance()
            .Build();

        var halfIntDescendants = registry.DescendantsOf<HalfIntegerMirrorClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F71MirrorSymmetryPi2Inheritance), halfIntDescendants);
    }

    [Fact]
    public void RegisterF71_AncestorsTransitivelyReachPolynomialFoundation()
    {
        var registry = BuildBaseRegistry()
            .RegisterF71MirrorSymmetryPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F71MirrorSymmetryPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Theory]
    [InlineData(3, 0, 1)]
    [InlineData(4, 0, 2)]
    [InlineData(5, 1, 2)]
    public void RegisterF71_MirrorPairAcrossRegistry(int N, int b, int expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF71MirrorSymmetryPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F71MirrorSymmetryPi2Inheritance>().MirrorPair(N, b));
    }

    [Theory]
    [InlineData(3, false)]   // odd
    [InlineData(4, true)]    // even
    [InlineData(5, false)]   // odd
    [InlineData(6, true)]    // even
    public void RegisterF71_HasCenterBondAcrossRegistry(int N, bool expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF71MirrorSymmetryPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F71MirrorSymmetryPi2Inheritance>().HasCenterBond(N));
    }

    [Fact]
    public void RegisterF71_LandsOnLadderAnchorAtSmallN()
    {
        var registry = BuildBaseRegistry()
            .RegisterF71MirrorSymmetryPi2Inheritance()
            .Build();

        var f = registry.Get<F71MirrorSymmetryPi2Inheritance>();
        // N=2,3 → 1 = a_1 (self-mirror); N=4,5 → 2 = a_0
        Assert.Equal(1, f.LadderIndexForIndependentComponentCount(2));
        Assert.Equal(1, f.LadderIndexForIndependentComponentCount(3));
        Assert.Equal(0, f.LadderIndexForIndependentComponentCount(4));
        Assert.Equal(0, f.LadderIndexForIndependentComponentCount(5));
        Assert.Null(f.LadderIndexForIndependentComponentCount(6));   // 3 not on ladder
    }
}
