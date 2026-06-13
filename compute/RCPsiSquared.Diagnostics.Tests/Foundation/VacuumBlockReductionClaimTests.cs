using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>VacuumBlockReductionClaim (the N=5 birth-canal boundary = the |1-exc&gt;&lt;vac| (0,1)
/// Liouville sector block): (a) the claim's own honest-scope metadata (Tier1Derived, N=5 scope,
/// the open V-Effect identity and absent aromaticity thesis named in the Summary); (b) the
/// bit-exact N=5 anchor the claim banks, recomputed live through its witness; (c) registry wiring
/// (typed parent AbsorptionTheoremClaim). Witness twin: <see cref="SectorReductionWitness"/>
/// (inspect --root reduction); Python twin: <c>simulations/birth_canal_vacuum_block_verifier.py</c>.</summary>
public class VacuumBlockReductionClaimTests
{
    // ------------------------------------------------------------------
    // (a) honest-scope metadata
    // ------------------------------------------------------------------

    [Fact]
    public void Claim_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, VacuumBlockReductionClaim.Shared.Tier);
    }

    [Fact]
    public void Claim_Summary_NamesN5ScopeAndTheOpenIdentity()
    {
        var s = VacuumBlockReductionClaim.Shared.Summary;
        Assert.Contains("N=5", s);
        Assert.Contains("(0,1)", s);
        // honest scope is load-bearing: the V-Effect self-pair identity stays OPEN, not claimed.
        Assert.Contains("OPEN", s, System.StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void Claim_CitesItsWitnessAndVerifier_BothDirections()
    {
        var a = VacuumBlockReductionClaim.Shared.Anchor;
        Assert.Contains("SectorReductionWitness", a);                       // C# witness named for discoverability
        Assert.Contains("birth_canal_vacuum_block_verifier.py", a);        // the verifier
    }

    // ------------------------------------------------------------------
    // (b) the N=5 anchor the claim banks (recomputed live via the witness)
    // ------------------------------------------------------------------

    [Fact]
    public void Claim_BanksTheBitExactN5CanalAnchor()
    {
        // The (0,1) block reproduces the full-4^N boundary at N=5 (canal profile), bit-exact;
        // the claim's banked rates are the SectorReductionWitness numbers, pinned here directly.
        var canal = new[] { 0.25, 1.5, 1.5, 1.5, 0.25 };
        Assert.Equal(1.2482918643729715,
            SectorReductionWitness.VacBlockSlowest(5, 1.5, canal, TopologyKind.Chain), 6);
        Assert.Equal(4.0 / 3.0,
            SectorReductionWitness.VacBlockSlowest(5, 1000.0, canal, TopologyKind.Chain), 5);
    }

    [Fact]
    public void Claim_FlatGammaBlindness_IsAnalyticAtEveryN()
    {
        // DERIVED: at uniform gamma, L_(1,0) = -iQ h - 2 gamma I, -iQh anti-Hermitian -> Re = -2 gamma,
        // Q-invariant. The claim asserts this analytic blindness; checked here past N=5.
        foreach (int n in new[] { 5, 6 })
        {
            var uni = System.Linq.Enumerable.Repeat(1.0, n).ToArray();
            Assert.Equal(2.0, SectorReductionWitness.VacBlockSlowest(n, 1.5, uni, TopologyKind.Chain), 9);
            Assert.Equal(2.0, SectorReductionWitness.VacBlockSlowest(n, 1000.0, uni, TopologyKind.Chain), 9);
        }
    }

    // ------------------------------------------------------------------
    // (c) registration
    // ------------------------------------------------------------------

    [Fact]
    public void BuildDefault_ContainsVacuumBlockReductionClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<VacuumBlockReductionClaim>());
    }

    [Fact]
    public void Claim_TierIsTier1Derived_FromRegistry()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.Equal(Tier.Tier1Derived, registry.Get<VacuumBlockReductionClaim>().Tier);
    }

    [Fact]
    public void Claim_Ancestors_ContainAbsorptionTheorem()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<VacuumBlockReductionClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
    }
}
