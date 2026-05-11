using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F2bXyChainSpectrumPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterF66PoleModesPi2Inheritance()
            .RegisterF65XxChainSpectrumPi2Inheritance();

    [Fact]
    public void RegisterF2b_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF2bXyChainSpectrumPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F2bXyChainSpectrumPi2Inheritance>());
    }

    [Fact]
    public void RegisterF2b_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF2bXyChainSpectrumPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F2bXyChainSpectrumPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF2b_AncestorsContainF65AndPi2DyadicLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF2bXyChainSpectrumPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F2bXyChainSpectrumPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F65XxChainSpectrumPi2Inheritance), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Theory]
    [InlineData(3, 1.0, 1, 1.4142135623730951)]
    [InlineData(3, 1.0, 2, 0.0)]
    public void RegisterF2b_EigenvalueAcrossRegistry(int N, double J, int k, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF2bXyChainSpectrumPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F2bXyChainSpectrumPi2Inheritance>().Eigenvalue(N, J, k), precision: 12);
    }

    [Fact]
    public void RegisterF2b_EigenvectorMatchesF65AcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF2bXyChainSpectrumPi2Inheritance()
            .Build();

        var f2b = registry.Get<F2bXyChainSpectrumPi2Inheritance>();
        for (int N = 2; N <= 5; N++)
            for (int k = 1; k <= N; k++)
                for (int site = 0; site < N; site++)
                    Assert.True(f2b.EigenvectorMatchesF65(N, k, site));
    }

    [Fact]
    public void RegisterF2b_WithoutF65_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterAbsorptionTheoremClaim()
                .RegisterF66PoleModesPi2Inheritance()
                // Missing: RegisterF65XxChainSpectrumPi2Inheritance
                .RegisterF2bXyChainSpectrumPi2Inheritance()
                .Build());
    }
}
