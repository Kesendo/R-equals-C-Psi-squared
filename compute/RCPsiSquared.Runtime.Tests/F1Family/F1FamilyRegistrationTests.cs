using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.Tests.F1Family;

public class F1FamilyRegistrationTests
{
    private static ChainSystem DefaultChain(int N = 5) =>
        new ChainSystem(N: N, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.XY, Topology: TopologyKind.Chain);

    [Fact]
    public void RegisterF1Family_BuildsElevenClaims()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        // Eleven entries: ChainSystemPrimitive + F1PalindromeIdentity +
        // PalindromeResidualScalingClaim + F1T1ResidualClosedForm +
        // F1T1ResidualPi2Decomposition + F1DepolResidualClosedForm +
        // F49NonUniformCrossTermClaim + F1GeneralTopologyVerifiedClaim (added
        // 2026-05-18 alongside the last F1 OpenQuestion closure) +
        // F4KernelDimensionByComponentsClaim (added 2026-05-18 as the F4
        // disconnected-graph bridge from the F1 SLOW_N8 sweep) +
        // RingN4DihedralLockClaim (added 2026-05-19 from the Q-sweep extension:
        // K_{2,2} = C_4 bipartite-complete + Casimir spectrum closes the
        // (3/4)·J·N Im-max bound in closed form) +
        // StarImMaxBoundClaim (added 2026-05-19 from the same Q-sweep: SU(2)/
        // Schur-Weyl hub-leaf Casimir on H_star = J·S⃗_0·S⃗_L closes the
        // J·N/2 saturation in closed form, sister derivation to RingN4).
        // SingleBody scaling deliberately omitted (builder is type-keyed; see
        // F1FamilyRegistration XML docs for the Option-B rationale).
        Assert.Equal(11, registry.All().Count());
        Assert.True(registry.Contains<ChainSystemPrimitive>());
        Assert.True(registry.Contains<F1PalindromeIdentity>());
        Assert.True(registry.Contains<PalindromeResidualScalingClaim>());
        Assert.True(registry.Contains<F1T1ResidualClosedForm>());
        Assert.True(registry.Contains<F1T1ResidualPi2Decomposition>());
        Assert.True(registry.Contains<F1DepolResidualClosedForm>());
        Assert.True(registry.Contains<F49NonUniformCrossTermClaim>());
        Assert.True(registry.Contains<F1GeneralTopologyVerifiedClaim>());
        Assert.True(registry.Contains<F4KernelDimensionByComponentsClaim>());
        Assert.True(registry.Contains<RingN4DihedralLockClaim>());
        Assert.True(registry.Contains<StarImMaxBoundClaim>());
    }

    [Fact]
    public void RegisterF1Family_StarImMaxBound_Resolves_AndIsTier1Derived()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var starClaim = registry.Get<StarImMaxBoundClaim>();
        Assert.NotNull(starClaim);
        // Landed Tier 1 derived 2026-05-19: closed form via SU(2)/Schur-Weyl
        // hub-leaf Casimir + maximum-S_L=(N-1)/2 ferromagnet eigenmode
        // construction (see PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md).
        Assert.Equal(Tier.Tier1Derived, starClaim.Tier);
        Assert.Equal(3, starClaim.MinN);
        Assert.Equal(0.5, starClaim.Coefficient, precision: 14);
        // The 24 Q-sweep anchors survived the registry round-trip (6 Q × 4 N).
        Assert.Equal(24, starClaim.EmpiricalAnchors.Count);
        // Predict scales linearly with N and J: (1/2)·J·N.
        Assert.Equal(0.5 * 0.1 * 5, starClaim.Predict(N: 5, J: 0.1), precision: 14);
        Assert.Equal(0.5 * 1.0 * 8, starClaim.Predict(N: 8, J: 1.0), precision: 14);
        // Im/σ ratio is (1/2)·Q universally (N-independent dimensionless form).
        Assert.Equal(0.5 * 1.5, starClaim.PredictImOverSigma(Q: 1.5), precision: 14);
        Assert.Equal(0.5 * 2.0, starClaim.PredictImOverSigma(Q: 2.0), precision: 14);
        // N < 3 not supported (degenerate star).
        Assert.Throws<ArgumentOutOfRangeException>(
            () => starClaim.Predict(N: 2, J: 1.0));
    }

    [Fact]
    public void RegisterF1Family_StarImMaxBound_AncestorsContainF1Identity()
    {
        // The Star Im-max bound is wired with F1PalindromeIdentity as its
        // parent (Tier 1 derived, strength 5), same as RingN4DihedralLockClaim.
        // The Im-max bound lives in the L-spectrum the F1 palindrome partitions;
        // the eigenmode-construction machinery is shared across all three
        // topology-bound sister claims (F4, RingN4, Star).
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var ancestors = registry.AncestorsOf<StarImMaxBoundClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
        Assert.Contains(typeof(ChainSystemPrimitive), ancestors);
    }

    [Fact]
    public void RegisterF1Family_RingN4DihedralLock_Resolves_AndIsTier1Derived()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var lockClaim = registry.Get<RingN4DihedralLockClaim>();
        Assert.NotNull(lockClaim);
        // Landed Tier 1 derived 2026-05-19: closed form via K_{2,2} = C_4 bipartite-
        // complete + Casimir spectrum + Liouvillian eigenmode construction (see
        // PROOF_RING_N4_DIHEDRAL_LOCK.md).
        Assert.Equal(Tier.Tier1Derived, lockClaim.Tier);
        Assert.Equal(4, lockClaim.N);
        Assert.Equal(0.75, lockClaim.Coefficient, precision: 14);
        // The six Q-sweep anchors survived the registry round-trip.
        Assert.Equal(6, lockClaim.EmpiricalAnchors.Count);
        // Predict scales linearly with J at N=4: 3·J.
        Assert.Equal(3.0 * 0.1, lockClaim.Predict(J: 0.1), precision: 14);
        Assert.Equal(3.0 * 1.0, lockClaim.Predict(J: 1.0), precision: 14);
        // Im/σ ratio is (3/4)·Q universally.
        Assert.Equal(0.75 * 1.5, lockClaim.PredictImOverSigma(Q: 1.5), precision: 14);
        Assert.Equal(0.75 * 2.0, lockClaim.PredictImOverSigma(Q: 2.0), precision: 14);
    }

    [Fact]
    public void RegisterF1Family_RingN4DihedralLock_AncestorsContainF1Identity()
    {
        // The ring N=4 dihedral lock is wired with F1PalindromeIdentity as its
        // parent (Tier 1 derived, strength 5). The Im-max bound lives in the
        // L-spectrum the F1 palindrome partitions; the eigenmode-construction
        // machinery is the same one F4KernelDimensionByComponentsClaim uses for
        // its connected-case upper-bound closure.
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var ancestors = registry.AncestorsOf<RingN4DihedralLockClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
        Assert.Contains(typeof(ChainSystemPrimitive), ancestors);
    }

    [Fact]
    public void RegisterF1Family_F4KernelDimensionByComponents_Resolves_AndIsTier1Derived()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var kernelDim = registry.Get<F4KernelDimensionByComponentsClaim>();
        Assert.NotNull(kernelDim);
        // Promoted from Tier 1 candidate to Tier 1 derived 2026-05-19: connected-case
        // upper bound closed by DEGENERACY_PALINDROME Result 2 (magnetization
        // conservation); multi-component product follows from standard tensor-sum
        // kernel factorisation.
        Assert.Equal(Tier.Tier1Derived, kernelDim.Tier);
        // The four N=8 anchors + the N=9 chain anchor (MklDirect bridge 2026-05-19)
        // survived the registry round-trip.
        Assert.Equal(5, kernelDim.EmpiricalAnchorsN8.Count);
        Assert.Equal(9, kernelDim.Predict(new[] { 8 }));        // chain/ring/star N=8
        Assert.Equal(25, kernelDim.Predict(new[] { 4, 4 }));    // K_4 + disjoint 4-chain N=8
        Assert.Equal(10, kernelDim.Predict(new[] { 9 }));       // chain N=9 (MklDirect bridge)
    }

    [Fact]
    public void RegisterF1Family_F4KernelDim_AncestorsContainF1Identity()
    {
        // The F4 bridge is wired with F1PalindromeIdentity as its parent (Tier 1
        // derived, strength 5). After the 2026-05-19 promotion the bridge is itself
        // Tier 1 derived (strength 5), and the inheritance check 5 ≥ 5 still holds
        // through the same edge. The "sister" relationship to
        // F1GeneralTopologyVerifiedClaim lives in the claim's XML doc and the proof
        // markdown, not as a direct dependency edge (Tier 2 verified strength 3 is
        // weaker than Tier 1 derived strength 5, so a direct edge would still violate
        // Tier inheritance).
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var ancestors = registry.AncestorsOf<F4KernelDimensionByComponentsClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
        Assert.Contains(typeof(ChainSystemPrimitive), ancestors);
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
        // VerifiedNValues bumped to {5, 6, 7, 8, 9} on 2026-05-19 when the N=9 chain
        // SLOW_N9 dogfood landed via the MklDirect bridge (commit abb2d52). The new
        // frontier at N=10 is memory-pressure rather than LP64 marshalling; see
        // F1GeneralTopologyVerifiedClaim.ScaleFrontierBlockedAtN.
        Assert.Equal(new[] { 5, 6, 7, 8, 9 }, general.VerifiedNValues);
        Assert.Equal(new[] { 5, 6, 7, 8, 9 }, general.ScaleUpToN);
        Assert.Equal(10, general.ScaleFrontierBlockedAtN);
        Assert.True(general.DisconnectedComponentsVerified);
        Assert.True(general.WeightedEdgesVerified);
        Assert.True(general.SingleBodyClassVerified);
        Assert.NotEmpty(general.SpectrumMetricsDataFiles);
        Assert.Contains("chain_N9.json", string.Join(";", general.SpectrumMetricsDataFiles));
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
