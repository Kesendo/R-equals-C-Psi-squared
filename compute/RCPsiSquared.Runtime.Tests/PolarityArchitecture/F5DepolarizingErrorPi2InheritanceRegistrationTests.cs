using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F5DepolarizingErrorPi2InheritanceRegistrationTests
{
    private static ChainSystem DefaultChain() =>
        new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.Heisenberg, Topology: TopologyKind.Chain);

    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror();

    [Fact]
    public void RegisterF5DepolarizingErrorPi2Inheritance_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF5DepolarizingErrorPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F5DepolarizingErrorPi2Inheritance>());
    }

    [Fact]
    public void RegisterF5DepolarizingErrorPi2Inheritance_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF5DepolarizingErrorPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F5DepolarizingErrorPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF5DepolarizingErrorPi2Inheritance_AncestorsContainAllThreeParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF5DepolarizingErrorPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F5DepolarizingErrorPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(Pi2OperatorSpaceMirrorClaim), ancestors);
    }

    [Fact]
    public void RegisterF5DepolarizingErrorPi2Inheritance_TwoOverThreeViaRegistry()
    {
        // Cross-registry verification: the Pi2-derived 2/3 ratio agrees with the
        // d/(d²-1) composition from the typed parents.
        var registry = BuildBaseRegistry()
            .RegisterF5DepolarizingErrorPi2Inheritance()
            .Build();

        var f = registry.Get<F5DepolarizingErrorPi2Inheritance>();
        Assert.Equal(2.0 / 3.0, f.TwoOverThree, precision: 12);
        Assert.Equal(2.0, f.DCoefficient, precision: 14);
        Assert.Equal(3.0, f.DSquaredMinusOne, precision: 14);
    }

    [Theory]
    [InlineData(3, 2.0 / 3.0)]
    [InlineData(4, 4.0 / 3.0)]
    [InlineData(5, 2.0)]
    public void RegisterF5DepolarizingErrorPi2Inheritance_LiveCoefficientAgreesWithClosedForm(int N, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF5DepolarizingErrorPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F5DepolarizingErrorPi2Inheritance>().LiveCoefficient(N), precision: 12);
    }

    [Fact]
    public void RegisterF5DepolarizingErrorPi2Inheritance_WithoutF1Family_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterF88PopcountCoherence()
                .RegisterF88StaticDyadicAnchor()
                .RegisterPi2OperatorSpaceMirror()
                // Missing: RegisterF1Family
                .RegisterF5DepolarizingErrorPi2Inheritance()
                .Build());
    }
}
