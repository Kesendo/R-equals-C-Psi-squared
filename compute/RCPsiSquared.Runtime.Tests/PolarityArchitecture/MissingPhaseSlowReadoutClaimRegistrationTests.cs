using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class MissingPhaseSlowReadoutClaimRegistrationTests
{
    private static readonly Lazy<RCPsiSquared.Runtime.ObjectManager.ClaimRegistry> Registry =
        new(KnowledgeRegistryFactory.BuildDefault);

    [Fact]
    public void RegistryOwnsTheReadoutWithExactlyItsThreePremises()
    {
        var registry = Registry.Value;
        var claim = registry.Get<MissingPhaseSlowReadoutClaim>();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        Assert.Same(registry.Get<MissingPhaseRelaxationScaleClaim>(), claim.RelaxationScale);
        Assert.Same(registry.Get<F70DeltaNSelectionRulePi2Inheritance>(), claim.SelectionRule);
        Assert.Same(registry.Get<ChiralKClaim>(), claim.Chirality);
        var parents = registry.EdgesInto<MissingPhaseSlowReadoutClaim>().Select(edge => edge.Parent).ToHashSet();
        Assert.True(parents.SetEquals(new[]
        {
            typeof(MissingPhaseRelaxationScaleClaim),
            typeof(F70DeltaNSelectionRulePi2Inheritance),
            typeof(ChiralKClaim)
        }));
        Assert.False(claim.SelectionRule.IsAllowedByDeltaNBound(2, 5));
    }

    [Fact]
    public void ReadoutInheritsGeometryAndLightButNotVacuumReduction()
    {
        var ancestors = Registry.Value.AncestorsOf<MissingPhaseSlowReadoutClaim>()
            .Select(claim => claim.GetType()).ToHashSet();
        Assert.Contains(typeof(JDefectLightMigrationClaim), ancestors);
        Assert.Contains(typeof(SeatCutBlindnessClaim), ancestors);
        Assert.Contains(typeof(JointPopcountSectors), ancestors);
        Assert.DoesNotContain(typeof(VacuumBlockReductionClaim), ancestors);
    }

    [Fact]
    public void ClaimLinksTheLiveEvidenceAndKeepsTheNonlinearReferenceFence()
    {
        var claim = Registry.Value.Get<MissingPhaseSlowReadoutClaim>();
        Assert.Contains("MissingPhaseSlowReadoutWitness.cs", claim.Anchor);
        Assert.Contains("PROOF_MISSING_PHASE_SLOW_READOUT.md", claim.Anchor);
        Assert.Contains("f(0)=-1/2", claim.PreparedSignal);
        Assert.Contains("d_out/d_2", claim.Scope);
        Assert.Contains("moving unitary reference", claim.Scope);
        Assert.Contains("No gamma-uniform", claim.Scope);
    }

    [Fact]
    public void ClaimRejectsMissingPremises()
    {
        var claim = Registry.Value.Get<MissingPhaseSlowReadoutClaim>();
        Assert.Throws<ArgumentNullException>(() => new MissingPhaseSlowReadoutClaim(null!, claim.SelectionRule, claim.Chirality));
        Assert.Throws<ArgumentNullException>(() => new MissingPhaseSlowReadoutClaim(claim.RelaxationScale, null!, claim.Chirality));
        Assert.Throws<ArgumentNullException>(() => new MissingPhaseSlowReadoutClaim(claim.RelaxationScale, claim.SelectionRule, null!));
    }
}
