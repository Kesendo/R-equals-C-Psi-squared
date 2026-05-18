using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F1DepolResidualClosedFormPi2InheritanceRegistrationTests
{
    private static ChainSystem DefaultChain() =>
        new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.Heisenberg, Topology: TopologyKind.Chain);

    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88bPopcountCoherence()
            .RegisterF88bStaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror();

    [Fact]
    public void RegisterF1DepolResidualPi2Inheritance_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF1DepolResidualClosedFormPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F1DepolResidualClosedFormPi2Inheritance>());
    }

    [Fact]
    public void RegisterF1DepolResidualPi2Inheritance_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF1DepolResidualClosedFormPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F1DepolResidualClosedFormPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF1DepolResidualPi2Inheritance_AncestorsContainAllThreeParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF1DepolResidualClosedFormPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F1DepolResidualClosedFormPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(Pi2OperatorSpaceMirrorClaim), ancestors);
    }

    [Fact]
    public void RegisterF1DepolResidualPi2Inheritance_CoefficientsViaRegistry()
    {
        // Cross-registry verification: the Pi2-derived coefficients agree bit-exact
        // with the parent closed form's constants. Drift on either side surfaces here.
        var registry = BuildBaseRegistry()
            .RegisterF1DepolResidualClosedFormPi2Inheritance()
            .Build();

        var f = registry.Get<F1DepolResidualClosedFormPi2Inheritance>();
        Assert.Equal(16.0 / 9.0, f.LocalCoefficient, precision: 14);
        Assert.Equal(16.0, f.CrossSiteCoefficient, precision: 14);
        Assert.Equal(4.0 / 3.0, f.PerPauliDepolarizingRate, precision: 14);
        Assert.Equal(4.0, f.DSquared, precision: 14);
        Assert.Equal(3.0, f.DSquaredMinusOne, precision: 14);
    }

    [Fact]
    public void RegisterF1DepolResidualPi2Inheritance_WithoutF1Family_Throws()
    {
        // Missing-parent guard for F1PalindromeIdentity: F1Family is the only
        // upstream registration that brings F1PalindromeIdentity into the builder.
        // Pin the specific Rule string ("MissingParent") to catch a regression
        // where the same exception type would surface for an unrelated reason
        // (e.g. Cycle, AnchorFileMissing, TierInheritance).
        var ex = Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterF88bPopcountCoherence()
                .RegisterF88bStaticDyadicAnchor()
                .RegisterPi2OperatorSpaceMirror()
                // Missing: RegisterF1Family
                .RegisterF1DepolResidualClosedFormPi2Inheritance()
                .Build());

        Assert.Equal("MissingParent", ex.Rule);
        Assert.Contains("F1PalindromeIdentity", ex.Message);
    }

    [Fact]
    public void RegisterF1DepolResidualPi2Inheritance_WithoutPi2Mirror_Throws()
    {
        // Missing-parent guard for Pi2OperatorSpaceMirrorClaim: needed for the
        // DSquaredMinusOne = 3 anchor that flows into (16/9, 16) (DSquared itself
        // now flows from Pi2DyadicLadderClaim.Term(-1), but DSquaredMinusOne still
        // anchors on the mirror's operator-space-pair semantic).
        var ex = Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF1Family(DefaultChain())
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterF88bPopcountCoherence()
                .RegisterF88bStaticDyadicAnchor()
                // Missing: RegisterPi2OperatorSpaceMirror
                .RegisterF1DepolResidualClosedFormPi2Inheritance()
                .Build());

        Assert.Equal("MissingParent", ex.Rule);
        Assert.Contains("Pi2OperatorSpaceMirrorClaim", ex.Message);
    }
}
