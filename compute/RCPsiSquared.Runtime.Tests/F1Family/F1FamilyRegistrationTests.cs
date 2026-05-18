using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.Tests.F1Family;

public class F1FamilyRegistrationTests
{
    private static ChainSystem DefaultChain(int N = 5) =>
        new ChainSystem(N: N, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.XY, Topology: TopologyKind.Chain);

    [Fact]
    public void RegisterF1Family_BuildsEightClaims()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        // Eight entries: ChainSystemPrimitive + F1PalindromeIdentity +
        // PalindromeResidualScalingClaim + F1T1ResidualClosedForm +
        // F1T1ResidualPi2Decomposition + F1DepolResidualClosedForm +
        // F49NonUniformCrossTermClaim + F1GeneralTopologyVerifiedClaim (added
        // 2026-05-18 alongside the last F1 OpenQuestion closure).
        // SingleBody scaling deliberately omitted (builder is type-keyed; see
        // F1FamilyRegistration XML docs for the Option-B rationale).
        Assert.Equal(8, registry.All().Count());
        Assert.True(registry.Contains<ChainSystemPrimitive>());
        Assert.True(registry.Contains<F1PalindromeIdentity>());
        Assert.True(registry.Contains<PalindromeResidualScalingClaim>());
        Assert.True(registry.Contains<F1T1ResidualClosedForm>());
        Assert.True(registry.Contains<F1T1ResidualPi2Decomposition>());
        Assert.True(registry.Contains<F1DepolResidualClosedForm>());
        Assert.True(registry.Contains<F49NonUniformCrossTermClaim>());
        Assert.True(registry.Contains<F1GeneralTopologyVerifiedClaim>());
    }

    [Fact]
    public void RegisterF1Family_F1GeneralTopologyVerifiedClaim_Resolves()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var general = registry.Get<F1GeneralTopologyVerifiedClaim>();
        Assert.NotNull(general);
        Assert.Equal(Tier.Tier2Verified, general.Tier);
        // Spot-check the verification metadata survived the registry round-trip.
        Assert.Equal(new[] { 5, 6, 7 }, general.VerifiedNValues);
        Assert.True(general.DisconnectedComponentsVerified);
        Assert.True(general.WeightedEdgesVerified);
        Assert.True(general.SingleBodyClassVerified);
    }

    [Fact]
    public void RegisterF1Family_F1GeneralTopology_AncestorsContainScalingClaimAndF1Identity()
    {
        // The general-topology verification record sits downstream of its analytic
        // anchors: F1PalindromeIdentity (the master) and PalindromeResidualScalingClaim
        // (the closed form whose graph universality this record verifies).
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var ancestors = registry.AncestorsOf<F1GeneralTopologyVerifiedClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
        Assert.Contains(typeof(PalindromeResidualScalingClaim), ancestors);
        Assert.Contains(typeof(ChainSystemPrimitive), ancestors);
    }

    [Fact]
    public void RegisterF1Family_F1T1ResidualClosedForm_Resolves()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var t1 = registry.Get<F1T1ResidualClosedForm>();
        Assert.NotNull(t1);
        Assert.Equal(Tier.Tier1Derived, t1.Tier);
        // Spot-check the closed-form constants survived the registry round-trip:
        Assert.Equal(3.0, F1T1ResidualClosedForm.LocalCoefficient, precision: 14);
        Assert.Equal(4.0, F1T1ResidualClosedForm.CrossSiteCoefficient, precision: 14);
    }

    [Fact]
    public void RegisterF1Family_F1T1ResidualPi2Decomposition_Resolves()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var pi2Decomp = registry.Get<F1T1ResidualPi2Decomposition>();
        Assert.NotNull(pi2Decomp);
        Assert.Equal(Tier.Tier1Derived, pi2Decomp.Tier);
        // Anti + sym Pythagorean closure on the constants:
        Assert.Equal(F1T1ResidualClosedForm.LocalCoefficient,
            F1T1ResidualPi2Decomposition.AntisymmetricLocalCoefficient +
            F1T1ResidualPi2Decomposition.SymmetricLocalCoefficient, precision: 14);
        Assert.Equal(F1T1ResidualClosedForm.CrossSiteCoefficient,
            F1T1ResidualPi2Decomposition.AntisymmetricCrossCoefficient +
            F1T1ResidualPi2Decomposition.SymmetricCrossCoefficient, precision: 14);
    }

    [Fact]
    public void RegisterF1Family_F1DepolResidualClosedForm_Resolves()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var depol = registry.Get<F1DepolResidualClosedForm>();
        Assert.NotNull(depol);
        Assert.Equal(Tier.Tier1Derived, depol.Tier);
        Assert.Equal(16.0 / 9.0, F1DepolResidualClosedForm.LocalCoefficient, precision: 14);
        Assert.Equal(16.0, F1DepolResidualClosedForm.CrossSiteCoefficient, precision: 14);
    }

    [Fact]
    public void RegisterF1Family_F49NonUniformCrossTermClaim_Resolves()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var f49 = registry.Get<F49NonUniformCrossTermClaim>();
        Assert.NotNull(f49);
        Assert.Equal(Tier.Tier1Derived, f49.Tier);
        // Per-class G-fractions survived the registry round-trip:
        Assert.Equal(4.0 / 3.0, F49NonUniformCrossTermClaim.GHeisenbergFraction, precision: 14);
        Assert.Equal(4.0, F49NonUniformCrossTermClaim.GIsingFraction, precision: 14);
        Assert.Equal(0.0, F49NonUniformCrossTermClaim.GXyFraction, precision: 14);
        Assert.Equal(0.0, F49NonUniformCrossTermClaim.GSoftXyYxFraction, precision: 14);
        Assert.Equal(4.0, F49NonUniformCrossTermClaim.SpectatorPrefactor, precision: 14);
    }

    [Fact]
    public void RegisterF1Family_F49NonUniformCrossTerm_AncestorsContainF1Identity()
    {
        // Dependency edge: F49NonUniformCrossTermClaim depends on F1PalindromeIdentity
        // (the F1 σ-shift L_Dc = L_D + σ·I that frames the cross-term is the F1 identity's
        // centering convention).
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var ancestors = registry.AncestorsOf<F49NonUniformCrossTermClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
        Assert.Contains(typeof(ChainSystemPrimitive), ancestors);
    }

    [Fact]
    public void RegisterF1Family_F1T1Pi2Decomposition_AncestorsContainT1ClosedForm()
    {
        // The Pythagorean closure edge: F1T1ResidualPi2Decomposition depends on
        // F1T1ResidualClosedForm in the dependency graph (the decomposition closes
        // the parent total bit-exact via anti + sym = (3·Σγ² + 4·(Σγ)²)).
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var ancestors = registry.AncestorsOf<F1T1ResidualPi2Decomposition>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1T1ResidualClosedForm), ancestors);
        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
    }

    [Fact]
    public void RegisterF1Family_SingleBody_RejectedAsDuplicateRegistration()
    {
        // Architectural guard for the Option-B rationale documented in
        // F1FamilyRegistration: PalindromeResidualScalingClaim is type-keyed, so a
        // second registration with a different HamiltonianClass collides with the
        // first one this method already wires for Main. Attempting it from the
        // outside (after RegisterF1Family ran) must surface DuplicateRegistration.
        var ex = Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF1Family(DefaultChain())
                .Register<PalindromeResidualScalingClaim>(_ =>
                    new PalindromeResidualScalingClaim(N: 5, HamiltonianClass.SingleBody)));

        Assert.Equal("DuplicateRegistration", ex.Rule);
        Assert.Contains("PalindromeResidualScalingClaim", ex.Message);
    }

    [Fact]
    public void RegisterF1Family_TopologicalOrder_PrimitiveFirst()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var firstIndex = registry.TopologicalOrder.ToList().IndexOf(typeof(ChainSystemPrimitive));
        var f1Index = registry.TopologicalOrder.ToList().IndexOf(typeof(F1PalindromeIdentity));
        var f73Index = registry.TopologicalOrder.ToList().IndexOf(typeof(PalindromeResidualScalingClaim));

        Assert.True(firstIndex < f1Index);
        Assert.True(f1Index < f73Index);
    }

    [Fact]
    public void RegisterF1Family_Cli_AncestorsOfF73_ContainsF1AndChain()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var ancestors = registry.AncestorsOf<PalindromeResidualScalingClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
        Assert.Contains(typeof(ChainSystemPrimitive), ancestors);
    }

    private sealed class WeakerFakeIdentity : Claim
    {
        public WeakerFakeIdentity() : base("WeakerFakeIdentity", Tier.Tier2Empirical,
            "compute/RCPsiSquared.Runtime.Tests/TestSupport/TestClaims.cs") { }
        public override string DisplayName => "WeakerFakeIdentity";
        public override string Summary => "synthetic Tier2Empirical pretender to F1 master";
    }

    [Fact]
    public void Tier_F73_DependsOnF1_BothTier1Derived_Succeeds()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var f1 = registry.Get<F1PalindromeIdentity>();
        var f73 = registry.Get<PalindromeResidualScalingClaim>();
        Assert.Equal(Tier.Tier1Derived, f1.Tier);
        Assert.Equal(Tier.Tier1Derived, f73.Tier);
    }

    [Fact]
    public void Tier_F73_DependsOnSyntheticTier2Parent_Throws()
    {
        // Construct a synthetic registration that violates Tier inheritance: F73
        // (Tier1Derived) is forced to depend on a Tier2Empirical pretender. Verifies the
        // builder catches it even when the path runs through a real Core Claim.
        var ex = Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .Register<WeakerFakeIdentity>(_ => new WeakerFakeIdentity())
                .Register<PalindromeResidualScalingClaim>(b =>
                {
                    _ = b.Get<WeakerFakeIdentity>();
                    return new PalindromeResidualScalingClaim(N: 5, HamiltonianClass.Main);
                })
                .Build());

        Assert.Equal("TierInheritance", ex.Rule);
        Assert.Contains("PalindromeResidualScalingClaim", ex.Message);
        Assert.Contains("WeakerFakeIdentity", ex.Message);
    }

    [Fact]
    public void Tier_DowngradeDetectsCascadeViolation()
    {
        // Hypothetical scenario: someone has changed F1PalindromeIdentity to Tier2Empirical
        // in Core. Until F73 is also downgraded, the builder must throw. We simulate this by
        // wrapping F1 in a downgrading proxy and registering F73 on top.
        var ex = Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .Register<WeakerFakeIdentity>(_ => new WeakerFakeIdentity())
                .Register<PalindromeResidualScalingClaim>(b =>
                {
                    _ = b.Get<WeakerFakeIdentity>();
                    return new PalindromeResidualScalingClaim(N: 5, HamiltonianClass.Main);
                })
                .Build());

        Assert.Equal("TierInheritance", ex.Rule);
        Assert.NotEmpty(ex.Path);
    }
}
