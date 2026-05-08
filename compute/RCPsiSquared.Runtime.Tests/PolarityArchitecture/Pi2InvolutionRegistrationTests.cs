using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class Pi2InvolutionRegistrationTests
{
    private static ChainSystem DefaultChain() =>
        new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.Heisenberg, Topology: TopologyKind.Chain);

    [Fact]
    public void RegisterPi2Involution_AddsClaim()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Involution()
            .Build();

        Assert.True(registry.Contains<Pi2InvolutionClaim>());
    }

    [Fact]
    public void RegisterPi2Involution_TierIsTier1Derived()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Involution()
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<Pi2InvolutionClaim>().Tier);
    }

    [Fact]
    public void RegisterPi2Involution_AncestorsContainF1Palindrome()
    {
        // Pi2Involution is the squaring consequence of F1; F1 must surface as ancestor.
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Involution()
            .Build();

        var ancestors = registry.AncestorsOf<Pi2InvolutionClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
    }

    [Fact]
    public void RegisterPi2Involution_WithoutF1Family_Throws()
    {
        // F1 must be registered first; without it the consequence-edge cannot be drawn.
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Involution()
                .Build());
    }
}
