using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class IbmBlockCpsiHardwareTableRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family();

    [Fact]
    public void RegisterIbmBlockCpsiHardwareTable_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterIbmBlockCpsiHardwareTable()
            .Build();

        Assert.True(registry.Contains<IbmBlockCpsiHardwareTable>());
    }

    [Fact]
    public void RegisterIbmBlockCpsiHardwareTable_TierIsTier2Verified()
    {
        var registry = BuildBaseRegistry()
            .RegisterIbmBlockCpsiHardwareTable()
            .Build();

        Assert.Equal(Tier.Tier2Verified, registry.Get<IbmBlockCpsiHardwareTable>().Tier);
    }

    [Fact]
    public void RegisterIbmBlockCpsiHardwareTable_AncestorsContainQuarter()
    {
        var registry = BuildBaseRegistry()
            .RegisterIbmBlockCpsiHardwareTable()
            .Build();

        var ancestors = registry.AncestorsOf<IbmBlockCpsiHardwareTable>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(QuarterAsBilinearMaxvalClaim), ancestors);
    }

    [Fact]
    public void RegisterIbmBlockCpsiHardwareTable_WithoutPi2Family_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                // Missing: RegisterPi2Family
                .RegisterIbmBlockCpsiHardwareTable()
                .Build());
    }
}
