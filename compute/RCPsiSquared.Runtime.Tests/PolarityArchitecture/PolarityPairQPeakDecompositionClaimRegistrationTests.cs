using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class PolarityPairQPeakDecompositionClaimRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family();

    [Fact]
    public void RegisterPolarityPairQPeakDecompositionClaim_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterPolarityPairQPeakDecompositionClaim()
            .Build();

        Assert.True(registry.Contains<PolarityPairQPeakDecompositionClaim>());
    }

    [Fact]
    public void RegisterPolarityPairQPeakDecompositionClaim_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterPolarityPairQPeakDecompositionClaim()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<PolarityPairQPeakDecompositionClaim>().Tier);
    }

    [Fact]
    public void RegisterPolarityPairQPeakDecompositionClaim_AncestorsContainAllThreeParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterPolarityPairQPeakDecompositionClaim()
            .Build();

        var ancestors = registry.AncestorsOf<PolarityPairQPeakDecompositionClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(HalfAsStructuralFixedPointClaim), ancestors);
        Assert.Contains(typeof(PolarityLayerOriginClaim), ancestors);
        Assert.Contains(typeof(QubitDimensionalAnchorClaim), ancestors);
    }

    [Fact]
    public void RegisterPolarityPairQPeakDecompositionClaim_WithoutPi2Family_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                // Missing: RegisterPi2Family
                .RegisterPolarityPairQPeakDecompositionClaim()
                .Build());
    }
}
