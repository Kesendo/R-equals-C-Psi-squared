using System.Linq;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class NivenRationalityRootClaimRegistrationTests
{
    [Fact]
    public void DefaultRegistry_ContainsNivenRationalityRootClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<NivenRationalityRootClaim>());
    }

    [Fact]
    public void NivenRationalityRootClaim_ResolvesFromDefaultRegistry_WithSharedTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var claim = registry.Get<NivenRationalityRootClaim>();
        Assert.NotNull(claim);

        // the registration must wire the SHARED registry instances, not fresh ones
        Assert.Same(registry.Get<TopologyBandEdgeClaim>(), claim.BandEdge);
        Assert.Same(registry.Get<F65XxChainSpectrumPi2Inheritance>(), claim.Rates);
    }

    [Fact]
    public void NivenRationalityRootClaim_HasBothSeFaceParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var ancestors = registry.AncestorsOf<NivenRationalityRootClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(TopologyBandEdgeClaim), ancestors);
        Assert.Contains(typeof(F65XxChainSpectrumPi2Inheritance), ancestors);
    }
}
