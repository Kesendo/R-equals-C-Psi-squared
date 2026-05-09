using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class DickeSuperpositionQuarterPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder();

    [Fact]
    public void RegisterDickeSuperpositionQuarter_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterDickeSuperpositionQuarterPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<DickeSuperpositionQuarterPi2Inheritance>());
    }

    [Fact]
    public void RegisterDickeSuperpositionQuarter_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterDickeSuperpositionQuarterPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<DickeSuperpositionQuarterPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterDickeSuperpositionQuarter_AncestorsContainBothPi2Anchors()
    {
        var registry = BuildBaseRegistry()
            .RegisterDickeSuperpositionQuarterPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<DickeSuperpositionQuarterPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(QuarterAsBilinearMaxvalClaim), ancestors);
        Assert.Contains(typeof(HalfAsStructuralFixedPointClaim), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Fact]
    public void RegisterDickeSuperpositionQuarter_AmGmIdentityHolds()
    {
        // Cross-registry verification: SectorBalance² = QuarterCeiling (= 1/4) exact.
        var registry = BuildBaseRegistry()
            .RegisterDickeSuperpositionQuarterPi2Inheritance()
            .Build();

        var d = registry.Get<DickeSuperpositionQuarterPi2Inheritance>();
        Assert.Equal(0.25, d.QuarterCeiling, precision: 14);
        Assert.Equal(0.5, d.SectorBalance, precision: 14);
        Assert.Equal(d.QuarterCeiling, d.SectorBalanceSquared, precision: 14);
    }

    [Theory]
    [InlineData(5, 2)]
    [InlineData(8, 3)]
    public void RegisterDickeSuperpositionQuarter_LiveBlockCpsiIsExactlyQuarter(int N, int n)
    {
        // Theorem 1 mechanism through the registry: M_block · 1/(4·M_block) = 1/4 exact.
        var registry = BuildBaseRegistry()
            .RegisterDickeSuperpositionQuarterPi2Inheritance()
            .Build();

        var d = registry.Get<DickeSuperpositionQuarterPi2Inheritance>();
        Assert.Equal(0.25, d.LiveBlockCpsiAtZero(N, n), precision: 12);
    }
}
