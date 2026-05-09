using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F39DetPiPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror();

    [Fact]
    public void RegisterF39DetPiPi2Inheritance_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF39DetPiPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F39DetPiPi2Inheritance>());
    }

    [Fact]
    public void RegisterF39DetPiPi2Inheritance_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF39DetPiPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F39DetPiPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF39DetPiPi2Inheritance_AncestorsContainBothPi2Anchors()
    {
        var registry = BuildBaseRegistry()
            .RegisterF39DetPiPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F39DetPiPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(Pi2OperatorSpaceMirrorClaim), ancestors);
    }

    [Theory]
    [InlineData(1, -1)]    // N=1: det = -1
    [InlineData(2, 1)]     // N≥2: det = +1
    [InlineData(3, 1)]
    [InlineData(4, 1)]
    public void RegisterF39DetPiPi2Inheritance_DetPiAcrossRegistry(int N, int expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF39DetPiPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F39DetPiPi2Inheritance>().DetPi(N));
    }

    [Fact]
    public void RegisterF39DetPiPi2Inheritance_PowerFactorMatchesLadder()
    {
        // Cross-registry verification: PowerNMinus1Factor(N) bit-exact equals
        // Pi2DyadicLadderClaim.Term(3-2N) for the same N.
        var registry = BuildBaseRegistry()
            .RegisterF39DetPiPi2Inheritance()
            .Build();

        var f39 = registry.Get<F39DetPiPi2Inheritance>();
        var ladder = registry.Get<Pi2DyadicLadderClaim>();

        for (int N = 1; N <= 5; N++)
            Assert.Equal(ladder.Term(f39.LadderIndexFor(N)), f39.PowerNMinus1Factor(N), precision: 12);
    }

    [Fact]
    public void RegisterF39DetPiPi2Inheritance_WithoutOperatorSpaceMirror_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                // Missing: F88* + Pi2OperatorSpaceMirror
                .RegisterF39DetPiPi2Inheritance()
                .Build());
    }
}
