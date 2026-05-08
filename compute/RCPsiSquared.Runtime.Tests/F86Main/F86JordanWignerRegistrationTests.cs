using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.F86Main;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.Tests.F86Main;

public class F86JordanWignerRegistrationTests
{
    [Fact]
    public void RegisterF86JordanWignerLight_BuildsFiveClaims()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86JordanWignerLight(N: 5, n: 1, gammaZero: 0.05)
            .Build();

        Assert.Equal(5, registry.All().Count());
        Assert.True(registry.Contains<XyJordanWignerModes>());
        Assert.True(registry.Contains<BondHdChannelWeights>());
        Assert.True(registry.Contains<JwBlockBasis>());
        Assert.True(registry.Contains<JwDispersionStructure>());
        Assert.True(registry.Contains<JwClusterDEigenstructure>());
    }

    [Fact]
    public void RegisterF86JordanWignerLight_AllPrimitivesAreTier1Derived()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86JordanWignerLight(N: 5, n: 1, gammaZero: 0.05)
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<XyJordanWignerModes>().Tier);
        Assert.Equal(Tier.Tier1Derived, registry.Get<BondHdChannelWeights>().Tier);
        Assert.Equal(Tier.Tier1Derived, registry.Get<JwBlockBasis>().Tier);
        Assert.Equal(Tier.Tier1Derived, registry.Get<JwDispersionStructure>().Tier);
        Assert.Equal(Tier.Tier1Derived, registry.Get<JwClusterDEigenstructure>().Tier);
    }

    [Fact]
    public void RegisterF86JordanWignerLight_JwBlockBasis_DescendsFromXyModes()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86JordanWignerLight(N: 5, n: 1, gammaZero: 0.05)
            .Build();

        var ancestors = registry.AncestorsOf<JwBlockBasis>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(XyJordanWignerModes), ancestors);
    }

    [Fact]
    public void RegisterF86JordanWignerLight_JwDispersionStructure_DescendsFromXyModes()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86JordanWignerLight(N: 5, n: 1, gammaZero: 0.05)
            .Build();

        var ancestors = registry.AncestorsOf<JwDispersionStructure>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(XyJordanWignerModes), ancestors);
    }

    [Fact]
    public void RegisterF86JordanWignerLight_JwClusterDEigenstructure_DescendsTransitivelyFromXyModes()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86JordanWignerLight(N: 5, n: 1, gammaZero: 0.05)
            .Build();

        var ancestors = registry.AncestorsOf<JwClusterDEigenstructure>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(JwBlockBasis), ancestors);
        Assert.Contains(typeof(JwDispersionStructure), ancestors);
        Assert.Contains(typeof(XyJordanWignerModes), ancestors);
    }

    [Fact]
    public void RegisterF86JordanWignerLight_BondHdChannelWeights_HasNoJwAncestors()
    {
        // Independent foundation: uses BlockLDecomposition.MhPerBond, not the JW modes.
        var registry = new ClaimRegistryBuilder()
            .RegisterF86JordanWignerLight(N: 5, n: 1, gammaZero: 0.05)
            .Build();

        var ancestors = registry.AncestorsOf<BondHdChannelWeights>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.DoesNotContain(typeof(XyJordanWignerModes), ancestors);
        Assert.DoesNotContain(typeof(JwBlockBasis), ancestors);
        Assert.DoesNotContain(typeof(JwDispersionStructure), ancestors);
    }

    [Fact]
    public void RegisterF86JordanWignerLight_AnchorsReferenceProofF86Qpeak()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86JordanWignerLight(N: 5, n: 1, gammaZero: 0.05)
            .Build();

        Assert.Contains("PROOF_F86_QPEAK.md", registry.Get<XyJordanWignerModes>().Anchor);
        Assert.Contains("PROOF_F86_QPEAK.md", registry.Get<JwBlockBasis>().Anchor);
        Assert.Contains("PROOF_F86_QPEAK.md", registry.Get<JwDispersionStructure>().Anchor);
        Assert.Contains("PROOF_F86_QPEAK.md", registry.Get<JwClusterDEigenstructure>().Anchor);
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void RegisterF86JordanWignerLight_AcrossN_BuildsCleanly(int N)
    {
        // Smoke across the JW-track's empirical N range. The build performs anchor-file
        // existence + tier-inheritance verification, so passing across N is the
        // architecture-level confirmation that the lightweight family is consistent.
        var registry = new ClaimRegistryBuilder()
            .RegisterF86JordanWignerLight(N, n: 1, gammaZero: 0.05)
            .Build();

        Assert.Equal(5, registry.All().Count());
    }

    [Fact]
    public void RegisterF86JordanWignerLight_NonC2Block_Throws()
    {
        // c=3 (n=2) is rejected by the c=2-only primitives downstream.
        Assert.Throws<ArgumentException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF86JordanWignerLight(N: 5, n: 2, gammaZero: 0.05)
                .Build());
    }
}
