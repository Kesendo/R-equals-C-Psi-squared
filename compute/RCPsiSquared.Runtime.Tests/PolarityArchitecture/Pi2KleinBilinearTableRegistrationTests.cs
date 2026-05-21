using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class Pi2KleinBilinearTableRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family();

    [Fact]
    public void RegisterPi2KleinBilinearTable_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterPi2KleinBilinearTable()
            .Build();

        Assert.True(registry.Contains<Pi2KleinBilinearTable>());
    }

    [Fact]
    public void RegisterPi2KleinBilinearTable_TierIsTier2Empirical()
    {
        var registry = BuildBaseRegistry()
            .RegisterPi2KleinBilinearTable()
            .Build();

        Assert.Equal(Tier.Tier2Empirical, registry.Get<Pi2KleinBilinearTable>().Tier);
    }

    [Fact]
    public void RegisterPi2KleinBilinearTable_AncestorsContainKleinFourCell()
    {
        var registry = BuildBaseRegistry()
            .RegisterPi2KleinBilinearTable()
            .Build();

        var ancestors = registry.AncestorsOf<Pi2KleinBilinearTable>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(KleinFourCellClaim), ancestors);
    }

    [Fact]
    public void RegisterPi2KleinBilinearTable_WithoutPi2Family_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                // Missing: RegisterPi2Family
                .RegisterPi2KleinBilinearTable()
                .Build());
    }
}
