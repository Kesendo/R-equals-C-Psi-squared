using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Schicht-1 wiring tests for the three F89 descendant claims:
/// F89AdditiveIdentity, F89PathKVacSeParseval, F89Path2Cardano.
/// All three share F89TopologyOrbitClosure as their parent edge.</summary>
public class F89DescendantClaimsRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterHalfIntegerMirror(N: 5)
            .RegisterF70DeltaNSelectionRulePi2Inheritance()
            .RegisterF72BlockDiagonalPurityPi2Inheritance()
            .RegisterF73SpatialSumPurityClosurePi2Inheritance()
            .RegisterF71MirrorSymmetryPi2Inheritance()
            .RegisterF89TopologyOrbitClosure();

    [Fact]
    public void RegisterF89AdditiveIdentity_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89AdditiveIdentityClaim()
            .Build();
        Assert.True(registry.Contains<F89AdditiveIdentityClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<F89AdditiveIdentityClaim>().Tier);
    }

    [Fact]
    public void RegisterF89AdditiveIdentity_AncestorsContainF89()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89AdditiveIdentityClaim()
            .Build();
        var ancestors = registry.AncestorsOf<F89AdditiveIdentityClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F89TopologyOrbitClosure), ancestors);
        // Transitive: should reach AbsorptionTheoremClaim and Pi2DyadicLadderClaim via F89
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Fact]
    public void RegisterF89PathKVacSeParseval_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89PathKVacSeParsevalClaim()
            .Build();
        Assert.True(registry.Contains<F89PathKVacSeParsevalClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<F89PathKVacSeParsevalClaim>().Tier);
    }

    [Fact]
    public void RegisterF89PathKVacSeParseval_AncestorsContainF89()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89PathKVacSeParsevalClaim()
            .Build();
        var ancestors = registry.AncestorsOf<F89PathKVacSeParsevalClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F89TopologyOrbitClosure), ancestors);
    }

    [Fact]
    public void RegisterF89Path2Cardano_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89Path2CardanoClaim()
            .Build();
        Assert.True(registry.Contains<F89Path2CardanoClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<F89Path2CardanoClaim>().Tier);
    }

    [Fact]
    public void RegisterF89Path2Cardano_AncestorsContainF89()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89Path2CardanoClaim()
            .Build();
        var ancestors = registry.AncestorsOf<F89Path2CardanoClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F89TopologyOrbitClosure), ancestors);
    }

    [Fact]
    public void RegisterAll_F89DescendantsAppearTogether()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89AdditiveIdentityClaim()
            .RegisterF89PathKVacSeParsevalClaim()
            .RegisterF89Path2CardanoClaim()
            .Build();
        var descendants = registry.DescendantsOf<F89TopologyOrbitClosure>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F89AdditiveIdentityClaim), descendants);
        Assert.Contains(typeof(F89PathKVacSeParsevalClaim), descendants);
        Assert.Contains(typeof(F89Path2CardanoClaim), descendants);
    }

    [Fact]
    public void RegisterF89AdditiveIdentity_BarePerSite_AccessibleAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89AdditiveIdentityClaim()
            .Build();
        // S_bare(t=0, N=7) = 6/49 (static method, but verify Claim is accessible)
        Assert.True(registry.Contains<F89AdditiveIdentityClaim>());
        Assert.Equal(6.0 / 49.0, F89AdditiveIdentityClaim.BarePerSite(7, 0.05, 0), precision: 14);
    }

    [Fact]
    public void RegisterF89Path2Cardano_VerifyAtStandardQ_AccessibleAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89Path2CardanoClaim()
            .Build();
        Assert.True(F89Path2CardanoClaim.VerifyAtStandardQ());
    }

    [Fact]
    public void RegisterF89PathKAtLockMechanism_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89PathKAtLockMechanismClaim()
            .Build();
        Assert.True(registry.Contains<F89PathKAtLockMechanismClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<F89PathKAtLockMechanismClaim>().Tier);
    }

    [Fact]
    public void RegisterF89UnifiedFaClosedForm_AddsClaim_AfterAtLock()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89PathKAtLockMechanismClaim()
            .RegisterF89UnifiedFaClosedFormClaim()
            .Build();
        Assert.True(registry.Contains<F89UnifiedFaClosedFormClaim>());
        // Should have BOTH F89 and F89PathKAtLockMechanism as ancestors
        var ancestors = registry.AncestorsOf<F89UnifiedFaClosedFormClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F89TopologyOrbitClosure), ancestors);
        Assert.Contains(typeof(F89PathKAtLockMechanismClaim), ancestors);
    }

    [Fact]
    public void RegisterF89Path3SeDeFactorisation_AddsClaim_AfterAtLock()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89PathKAtLockMechanismClaim()
            .RegisterF89Path3SeDeFactorisationClaim()
            .Build();
        Assert.True(registry.Contains<F89Path3SeDeFactorisationClaim>());
        var ancestors = registry.AncestorsOf<F89Path3SeDeFactorisationClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F89TopologyOrbitClosure), ancestors);
        Assert.Contains(typeof(F89PathKAtLockMechanismClaim), ancestors);
    }

    [Fact]
    public void RegisterF89Path3OcticEp_AddsClaim_AfterAtLock()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89PathKAtLockMechanismClaim()
            .RegisterF89Path3OcticEpClaim()
            .Build();
        Assert.True(registry.Contains<F89Path3OcticEpClaim>());
        var ancestors = registry.AncestorsOf<F89Path3OcticEpClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F89PathKAtLockMechanismClaim), ancestors);
    }

    [Fact]
    public void RegisterF89Path3OcticGalois_AddsClaim_AfterAtLock()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89PathKAtLockMechanismClaim()
            .RegisterF89Path3OcticGaloisClaim()
            .Build();
        Assert.True(registry.Contains<F89Path3OcticGaloisClaim>());
        var ancestors = registry.AncestorsOf<F89Path3OcticGaloisClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F89PathKAtLockMechanismClaim), ancestors);
    }

    [Fact]
    public void RegisterF89PathKHbMixedDegrees_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89PathKHbMixedDegreesClaim()
            .Build();
        Assert.True(registry.Contains<F89PathKHbMixedDegreesClaim>());
        var ancestors = registry.AncestorsOf<F89PathKHbMixedDegreesClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F89TopologyOrbitClosure), ancestors);
    }

    [Fact]
    public void RegisterAll_NewSixClaims_AppearTogetherUnderF89()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89AdditiveIdentityClaim()
            .RegisterF89PathKVacSeParsevalClaim()
            .RegisterF89Path2CardanoClaim()
            .RegisterF89PathKAtLockMechanismClaim()
            .RegisterF89UnifiedFaClosedFormClaim()
            .RegisterF89Path3SeDeFactorisationClaim()
            .RegisterF89Path3OcticEpClaim()
            .RegisterF89Path3OcticGaloisClaim()
            .RegisterF89PathKHbMixedDegreesClaim()
            .Build();
        var descendants = registry.DescendantsOf<F89TopologyOrbitClosure>()
            .Select(c => c.GetType()).ToHashSet();
        // Original three
        Assert.Contains(typeof(F89AdditiveIdentityClaim), descendants);
        Assert.Contains(typeof(F89PathKVacSeParsevalClaim), descendants);
        Assert.Contains(typeof(F89Path2CardanoClaim), descendants);
        // New six
        Assert.Contains(typeof(F89PathKAtLockMechanismClaim), descendants);
        Assert.Contains(typeof(F89UnifiedFaClosedFormClaim), descendants);
        Assert.Contains(typeof(F89Path3SeDeFactorisationClaim), descendants);
        Assert.Contains(typeof(F89Path3OcticEpClaim), descendants);
        Assert.Contains(typeof(F89Path3OcticGaloisClaim), descendants);
        Assert.Contains(typeof(F89PathKHbMixedDegreesClaim), descendants);
    }
}
