using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class MissingPhaseRelaxationScaleClaimRegistrationTests
{
    private static readonly Lazy<RCPsiSquared.Runtime.ObjectManager.ClaimRegistry> Registry =
        new(KnowledgeRegistryFactory.BuildDefault);

    [Fact]
    public void DefaultRegistry_ContainsTheTier1DerivedClaim()
    {
        var registry = Registry.Value;

        Assert.True(registry.Contains<MissingPhaseRelaxationScaleClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<MissingPhaseRelaxationScaleClaim>().Tier);
    }

    [Fact]
    public void EdgesIntoClaim_AreExactlyTheThreeDirectPremises()
    {
        var parents = Registry.Value.EdgesInto<MissingPhaseRelaxationScaleClaim>()
            .Select(edge => edge.Parent)
            .ToHashSet();

        var expected = new HashSet<Type>
        {
            typeof(JDefectLightMigrationClaim),
            typeof(SeatCutBlindnessClaim),
            typeof(JointPopcountSectors)
        };
        Assert.True(parents.SetEquals(expected),
            $"expected [{string.Join(", ", expected.Select(type => type.Name))}], " +
            $"got [{string.Join(", ", parents.Select(type => type.Name))}]");
    }

    [Fact]
    public void RegisteredClaim_HoldsTheSameInstancesAsItsEdges()
    {
        var registry = Registry.Value;
        var claim = registry.Get<MissingPhaseRelaxationScaleClaim>();

        Assert.Same(registry.Get<JDefectLightMigrationClaim>(), claim.LightMigration);
        Assert.Same(registry.Get<SeatCutBlindnessClaim>(), claim.SeatBlindness);
        Assert.Same(registry.Get<JointPopcountSectors>(), claim.JointSectors);
    }

    [Fact]
    public void GeneralAbsorptionAndPalindromeAreAncestorsButVacuumReductionIsNot()
    {
        var registry = Registry.Value;
        var ancestors = registry.AncestorsOf<MissingPhaseRelaxationScaleClaim>()
            .Select(claim => claim.GetType())
            .ToHashSet();
        var mechanismParents = registry.EdgesInto<JDefectLightMigrationClaim>()
            .Select(edge => edge.Parent)
            .ToHashSet();

        Assert.Contains(typeof(JDefectLightMigrationClaim), ancestors);
        Assert.Contains(typeof(SeatCutBlindnessClaim), ancestors);
        Assert.Contains(typeof(JointPopcountSectors), ancestors);
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
        Assert.Contains(typeof(AbsorptionTheoremClaim), mechanismParents);
        Assert.Contains(typeof(F1PalindromeIdentity), mechanismParents);
        Assert.DoesNotContain(typeof(VacuumBlockReductionClaim), ancestors);
    }
}
