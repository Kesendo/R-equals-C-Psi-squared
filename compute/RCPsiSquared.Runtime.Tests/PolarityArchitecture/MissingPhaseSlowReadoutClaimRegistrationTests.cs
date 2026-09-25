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
        Assert.Contains($"{nameof(MissingPhaseSlowReadoutWitness)}.cs", claim.Anchor);
        Assert.Contains("PROOF_MISSING_PHASE_SLOW_READOUT.md", claim.Anchor);
        // The values the claim states are the ones the live witness computes.
        var reading = new MissingPhaseSlowReadoutWitness().Reading;
        Assert.Equal(7, reading.NullPairsAtUniformPoint.Count());   // the ties below are not vacuous
        Assert.Equal(3, reading.NullPairsAtFirstOrder.Count());
        Assert.Contains($"f(0)={reading.EndZzResidue}", claim.PreparedSignal);
        Assert.All(reading.NullPairsAtUniformPoint, pair => Assert.Contains($"({pair.A},{pair.B})", claim.ReadoutZeros));
        Assert.Contains(string.Join(", ", reading.NullPairsAtFirstOrder.Select(pair => $"({pair.A},{pair.B})")),
            claim.ReadoutZeros);
        // The scope fences of a Tier-1 claim.
        Assert.Contains("d_out/d_2", claim.Scope);
        Assert.Contains("moving unitary reference", claim.Scope);
        Assert.Contains("No gamma-uniform", claim.Scope);
        Assert.Contains("claimed at epsilon=0 only", claim.Scope);
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
