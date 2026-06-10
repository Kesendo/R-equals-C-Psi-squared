using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F86Main;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.F86Main;

public class F86JordanWignerRegistrationTests
{
    /// <summary>The JW light family + its ChiralKClaim parent (typed chirality edge,
    /// 2026-06-10: the dispersion is ChiralKClaim's BDI spectrum inversion
    /// ε_{N+1−k} = −ε_k, so ChiralKClaim must be registered in the same builder).</summary>
    private static ClaimRegistry BuildJwRegistry(int N = 5, int n = 1) =>
        new ClaimRegistryBuilder()
            .RegisterChiralK()
            .RegisterF86JordanWignerLight(N, n, gammaZero: 0.05)
            .Build();

    [Fact]
    public void RegisterF86JordanWignerLight_BuildsFiveJwClaims_PlusChiralKParent()
    {
        var registry = BuildJwRegistry();

        Assert.Equal(6, registry.All().Count());
        Assert.True(registry.Contains<XyJordanWignerModes>());
        Assert.True(registry.Contains<BondHdChannelWeights>());
        Assert.True(registry.Contains<JwBlockBasis>());
        Assert.True(registry.Contains<JwDispersionStructure>());
        Assert.True(registry.Contains<JwClusterDEigenstructure>());
        Assert.True(registry.Contains<ChiralKClaim>());
    }

    [Fact]
    public void RegisterF86JordanWignerLight_AllPrimitivesAreTier1Derived()
    {
        var registry = BuildJwRegistry();

        Assert.Equal(Tier.Tier1Derived, registry.Get<XyJordanWignerModes>().Tier);
        Assert.Equal(Tier.Tier1Derived, registry.Get<BondHdChannelWeights>().Tier);
        Assert.Equal(Tier.Tier1Derived, registry.Get<JwBlockBasis>().Tier);
        Assert.Equal(Tier.Tier1Derived, registry.Get<JwDispersionStructure>().Tier);
        Assert.Equal(Tier.Tier1Derived, registry.Get<JwClusterDEigenstructure>().Tier);
        // The chirality parent is at least as strong, so the typed edge keeps every tier
        // unchanged (TierInheritance invariant).
        Assert.Equal(Tier.Tier1Derived, registry.Get<ChiralKClaim>().Tier);
    }

    [Fact]
    public void RegisterF86JordanWignerLight_XyModes_DescendFromChiralK()
    {
        // The 2026-06-10 chirality edge: ε_{N+1−k} = −ε_k is ChiralKClaim's BDI spectrum
        // inversion, so the dispersion owner descends from the chiral K root.
        var registry = BuildJwRegistry();

        var ancestors = registry.AncestorsOf<XyJordanWignerModes>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(ChiralKClaim), ancestors);
    }

    [Fact]
    public void RegisterF86JordanWignerLight_ClusterTables_DescendTransitivelyFromChiralK()
    {
        // The cluster-difference tables inherit the ± pairing through XyJordanWignerModes.
        var registry = BuildJwRegistry();

        var dispersionAncestors = registry.AncestorsOf<JwDispersionStructure>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(ChiralKClaim), dispersionAncestors);

        var eigenAncestors = registry.AncestorsOf<JwClusterDEigenstructure>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(ChiralKClaim), eigenAncestors);
    }

    [Fact]
    public void RegisterF86JordanWignerLight_JwBlockBasis_DescendsFromXyModes()
    {
        var registry = BuildJwRegistry();

        var ancestors = registry.AncestorsOf<JwBlockBasis>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(XyJordanWignerModes), ancestors);
    }

    [Fact]
    public void RegisterF86JordanWignerLight_JwDispersionStructure_DescendsFromXyModes()
    {
        var registry = BuildJwRegistry();

        var ancestors = registry.AncestorsOf<JwDispersionStructure>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(XyJordanWignerModes), ancestors);
    }

    [Fact]
    public void RegisterF86JordanWignerLight_JwClusterDEigenstructure_DescendsTransitivelyFromXyModes()
    {
        var registry = BuildJwRegistry();

        var ancestors = registry.AncestorsOf<JwClusterDEigenstructure>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(JwBlockBasis), ancestors);
        Assert.Contains(typeof(JwDispersionStructure), ancestors);
        Assert.Contains(typeof(XyJordanWignerModes), ancestors);
    }

    [Fact]
    public void RegisterF86JordanWignerLight_BondHdChannelWeights_HasNoJwAncestors()
    {
        // Independent foundation: uses BlockLDecomposition.MhPerBond, not the JW modes,
        // so it inherits neither the JW family nor the chirality edge.
        var registry = BuildJwRegistry();

        var ancestors = registry.AncestorsOf<BondHdChannelWeights>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.DoesNotContain(typeof(XyJordanWignerModes), ancestors);
        Assert.DoesNotContain(typeof(JwBlockBasis), ancestors);
        Assert.DoesNotContain(typeof(JwDispersionStructure), ancestors);
        Assert.DoesNotContain(typeof(ChiralKClaim), ancestors);
    }

    [Fact]
    public void RegisterF86JordanWignerLight_AnchorsReferenceProofF86Qpeak()
    {
        var registry = BuildJwRegistry();

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
        var registry = BuildJwRegistry(N);

        Assert.Equal(6, registry.All().Count());
    }

    [Fact]
    public void RegisterF86JordanWignerLight_NonC2Block_Throws()
    {
        // c=3 (n=2) is rejected by the c=2-only primitives downstream.
        Assert.Throws<ArgumentException>(() => BuildJwRegistry(N: 5, n: 2));
    }
}
