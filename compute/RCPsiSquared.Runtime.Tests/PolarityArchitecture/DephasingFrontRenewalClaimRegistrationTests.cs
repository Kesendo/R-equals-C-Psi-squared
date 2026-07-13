using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>The dephasing-front renewal representation, typed as <see cref="DephasingFrontRenewalClaim"/>
/// (Tier1Derived; the exact solution of the watched walk). Wiring guard: the claim resolves from the default
/// registry with its two SHARED Tier1Derived typed parents (AbsorptionTheoremClaim = the sector rate Γ = 4γ,
/// F2bXyChainSpectrumPi2Inheritance = the clean propagator / band), and its live battery all-passes.</summary>
public class DephasingFrontRenewalClaimRegistrationTests
{
    [Fact]
    public void DefaultRegistry_ContainsDephasingFrontRenewalClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<DephasingFrontRenewalClaim>());
    }

    [Fact]
    public void Claim_ResolvesWithSharedTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<DephasingFrontRenewalClaim>();
        Assert.NotNull(claim);
        Assert.Same(registry.Get<AbsorptionTheoremClaim>(), claim.RateLaw);
        Assert.Same(registry.Get<F2bXyChainSpectrumPi2Inheritance>(), claim.Band);
    }

    [Fact]
    public void Claim_HasBothTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<DephasingFrontRenewalClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.Contains(typeof(F2bXyChainSpectrumPi2Inheritance), ancestors);
    }

    [Fact]
    public void Claim_BatteryAllPass()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<DephasingFrontRenewalClaim>();
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }
}
