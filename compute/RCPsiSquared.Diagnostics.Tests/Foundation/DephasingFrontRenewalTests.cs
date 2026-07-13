using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>DephasingFrontRenewalClaim / DephasingFrontRenewalWitness (the exact renewal
/// representation of the watched walk, docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md).
/// (a) the witness's live battery all-passes (renewal-vs-RK4, probability conservation, the
/// j=0 coherent-front Bessel identity, the Γ=0 clean-wave limit, the Haken-Strobl plateau, the
/// I₁ Airy constant); (b) the claim carries the same battery; (c) registry wiring: the claim
/// resolves from the default registry, Tier1Derived, with its two typed parents (the Absorption
/// Theorem rate Γ=4γ + the F2b clean band/propagator).</summary>
public class DephasingFrontRenewalTests
{
    // ------------------------------------------------------------------
    // (a) the witness battery (recomputed live at inspect time)
    // ------------------------------------------------------------------

    [Fact]
    public void Witness_BatteryAllPass()
    {
        var w = new DephasingFrontRenewalWitness();
        Assert.NotEmpty(w.Cases);
        Assert.Equal(6, w.Cases.Count);
        foreach (var c in w.Cases)
            Assert.True(c.Passes, $"battery case '{c.Name}': expected {c.Expected}, got {c.Actual} ({c.Detail})");
    }

    [Fact]
    public void Witness_Summary_NamesTheRenewalRepresentation()
    {
        var w = new DephasingFrontRenewalWitness();
        Assert.Contains("renewal", w.Summary, StringComparison.OrdinalIgnoreCase);
    }

    // ------------------------------------------------------------------
    // (b) the claim's battery (delegates to the witness)
    // ------------------------------------------------------------------

    [Fact]
    public void Claim_BatteryAllPass()
    {
        var claim = KnowledgeRegistryFactory.BuildDefault().Get<DephasingFrontRenewalClaim>();
        Assert.NotEmpty(claim.Cases);
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }

    // ------------------------------------------------------------------
    // (c) registration + tier + typed parents
    // ------------------------------------------------------------------

    [Fact]
    public void BuildDefault_ContainsDephasingFrontRenewalClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<DephasingFrontRenewalClaim>());
    }

    [Fact]
    public void Claim_TierIsTier1Derived()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.Equal(Tier.Tier1Derived, registry.Get<DephasingFrontRenewalClaim>().Tier);
    }

    [Fact]
    public void Claim_ResolvesWithSharedTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<DephasingFrontRenewalClaim>();
        Assert.Same(registry.Get<AbsorptionTheoremClaim>(), claim.RateLaw);
        Assert.Same(registry.Get<F2bXyChainSpectrumPi2Inheritance>(), claim.Band);
    }

    [Fact]
    public void Claim_Ancestors_ContainBothTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<DephasingFrontRenewalClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.Contains(typeof(F2bXyChainSpectrumPi2Inheritance), ancestors);
    }
}
