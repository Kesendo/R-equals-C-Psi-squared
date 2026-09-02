using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.F1Family;

/// <summary>Wiring tests for <see cref="CrackedRingExactCurveClaim"/> (F160), whose two typed parents are
/// <see cref="TopologyBandEdgeClaim"/> (the ring row and the Scope fence that is this curve's Perron root) and
/// <see cref="F2bXyChainSpectrumPi2Inheritance"/> (the u = 0 end). Built on the default registry, since the
/// two parents come from two different families.</summary>
public class CrackedRingExactCurveClaimRegistrationTests
{
    [Fact]
    public void DefaultRegistry_HoldsTheClaim_Tier1Derived()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<CrackedRingExactCurveClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<CrackedRingExactCurveClaim>().Tier);
    }

    [Fact]
    public void BothParents_AreAncestors_AndTheSameInstancesAsTheTypedEdges()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<CrackedRingExactCurveClaim>();
        var ancestors = registry.AncestorsOf<CrackedRingExactCurveClaim>().Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(TopologyBandEdgeClaim), ancestors);
        Assert.Contains(typeof(F2bXyChainSpectrumPi2Inheritance), ancestors);
        Assert.Same(registry.Get<TopologyBandEdgeClaim>(), claim.BandEdge);
        Assert.Same(registry.Get<F2bXyChainSpectrumPi2Inheritance>(), claim.ChainEnd);
    }

    [Fact]
    public void TheStatement_CarriesTheParityOfTheCount_AndTheKnobsSign()
    {
        // The first draft of the departure count forgot the parity of N (CAUGHT_ERRORS 2026-08-31); the two
        // committed conventions of the knob collide (the crack weakens, the walk-time step strengthens), so the
        // statement is in u. Neither may be lost in a later edit.
        var name = KnowledgeRegistryFactory.BuildDefault().Get<CrackedRingExactCurveClaim>().Name;
        Assert.Contains("for even N", name);
        Assert.Contains("for odd N", name);
        Assert.Contains("(N+1)/(N-1)", name);
        Assert.Contains("u = J'/J", name);
    }
}
