using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F49Pi2InheritanceRegistrationTests
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
    public void RegisterF49Pi2Inheritance_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF49Pi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F49Pi2Inheritance>());
    }

    [Fact]
    public void RegisterF49Pi2Inheritance_AncestorsContainAllThreeParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF49Pi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F49Pi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PalindromeResidualScalingClaim), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(Pi2OperatorSpaceMirrorClaim), ancestors);
    }

    [Theory]
    [InlineData(3, 4.0)]
    [InlineData(4, 16.0)]
    [InlineData(5, 64.0)]
    [InlineData(8, 4096.0)]
    public void RegisterF49Pi2Inheritance_PowerFactor_AgreesAcrossRegistry(int chainN, double expected)
    {
        // Cross-registry verification: F49 PowerFactor(chainN) bit-exact equals the
        // operator-space pinned dimension for (chainN−2) qubits in the mirror table.
        var registry = BuildBaseRegistry()
            .RegisterF49Pi2Inheritance()
            .Build();

        var f = registry.Get<F49Pi2Inheritance>();
        Assert.Equal(expected, f.PowerFactor(chainN), precision: 12);
        Assert.Equal(expected, f.MirrorPinnedPowerFactor(chainN), precision: 12);
    }

    [Fact]
    public void RegisterF49Pi2Inheritance_WithoutOperatorSpaceMirror_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF1Family(DefaultChain())
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                // Missing: RegisterF88PopcountCoherence + RegisterF88StaticDyadicAnchor + RegisterPi2OperatorSpaceMirror
                .RegisterF49Pi2Inheritance()
                .Build());
    }
}
