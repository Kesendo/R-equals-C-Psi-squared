using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Diagnostics.Tests.F87;

public class F87Pi2InheritanceRegistrationTests
{
    private static ChainSystem DefaultChain() =>
        new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.Heisenberg, Topology: TopologyKind.Chain);

    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF1Pi2Inheritance()
            .RegisterF87Family();

    [Fact]
    public void RegisterF87Pi2Inheritance_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF87Pi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F87Pi2Inheritance>());
    }

    [Fact]
    public void RegisterF87Pi2Inheritance_AncestorsContainAllThreeParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF87Pi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F87Pi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F87TrichotomyClassification), ancestors);
        Assert.Contains(typeof(F1Pi2Inheritance), ancestors);
        Assert.Contains(typeof(KleinFourCellClaim), ancestors);
    }

    [Fact]
    public void RegisterF87Pi2Inheritance_TransitivelyInheritedTwoFactorIsTwo()
    {
        // The "2" in F87's discriminator M = F1's residual; via F1Pi2Inheritance the "2"
        // is a_0 on the Pi2 ladder. Cross-registry verification.
        var registry = BuildBaseRegistry()
            .RegisterF87Pi2Inheritance()
            .Build();

        Assert.Equal(2.0, registry.Get<F87Pi2Inheritance>().TransitivelyInheritedTwoFactor, precision: 14);
    }

    [Fact]
    public void RegisterF87Pi2Inheritance_WithoutF1Inheritance_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF1Family(DefaultChain())
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                // Missing: RegisterF1Pi2Inheritance
                .RegisterF87Family()
                .RegisterF87Pi2Inheritance()
                .Build());
    }
}
