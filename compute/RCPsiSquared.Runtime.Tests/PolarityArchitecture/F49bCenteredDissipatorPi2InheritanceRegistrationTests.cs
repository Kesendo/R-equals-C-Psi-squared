using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F49bCenteredDissipatorPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror();

    [Fact]
    public void RegisterF49bCenteredDissipatorPi2Inheritance_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF49bCenteredDissipatorPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F49bCenteredDissipatorPi2Inheritance>());
    }

    [Fact]
    public void RegisterF49bCenteredDissipatorPi2Inheritance_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF49bCenteredDissipatorPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F49bCenteredDissipatorPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF49bCenteredDissipatorPi2Inheritance_AncestorsContainBothPi2Anchors()
    {
        var registry = BuildBaseRegistry()
            .RegisterF49bCenteredDissipatorPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F49bCenteredDissipatorPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(Pi2OperatorSpaceMirrorClaim), ancestors);
    }

    [Theory]
    [InlineData(1, 4.0)]
    [InlineData(2, 16.0)]
    [InlineData(3, 64.0)]
    [InlineData(6, 4096.0)]
    public void RegisterF49bCenteredDissipatorPi2Inheritance_FourPowerNAcrossRegistry(int N, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF49bCenteredDissipatorPi2Inheritance()
            .Build();

        var f = registry.Get<F49bCenteredDissipatorPi2Inheritance>();
        Assert.Equal(expected, f.FourPowerNFactor(N), precision: 12);
        Assert.Equal(expected, f.MirrorPinnedFourPowerN(N), precision: 12);
    }

    [Fact]
    public void RegisterF49bCenteredDissipatorPi2Inheritance_LadderAndMirrorAgree()
    {
        // Direct cross-check: F49b's 4^N from the ladder agrees bit-exactly with
        // Pi2OperatorSpaceMirror's pinned d² for N qubits, across N=1..6.
        var registry = BuildBaseRegistry()
            .RegisterF49bCenteredDissipatorPi2Inheritance()
            .Build();

        var f = registry.Get<F49bCenteredDissipatorPi2Inheritance>();
        for (int N = 1; N <= 6; N++)
            Assert.Equal(f.FourPowerNFactor(N), f.MirrorPinnedFourPowerN(N), precision: 12);
    }

    [Fact]
    public void RegisterF49bCenteredDissipatorPi2Inheritance_WithoutOperatorSpaceMirror_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                // Missing: F88* + Pi2OperatorSpaceMirror
                .RegisterF49bCenteredDissipatorPi2Inheritance()
                .Build());
    }
}
