using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Orchestration.Sweep;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Orchestration.Tests.Sweep;

public class SweepCoordinatorTests
{
    private static ChainSystem DefaultChain(int N = 5) =>
        new ChainSystem(N: N, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.XY, Topology: TopologyKind.Chain);

    [Fact]
    public void Sweep_F73Frobenius_30Points_AllMatched()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var coord = new SweepCoordinator(registry);
        var result = coord.Sweep(new SweepDimension.F73Frobenius(
            NValues: new[] { 2, 3, 4, 5, 6 },
            HClasses: new[] { HamiltonianClass.Main, HamiltonianClass.SingleBody },
            ChainOnly: new[] { true, true, true })); // three γ surrogates; F73 closed form is γ-independent

        Assert.Equal(30, result.Total);
        Assert.Equal(30, result.Matched);
        Assert.Equal(0, result.Mismatched);
    }

    [Fact]
    public void Sweep_Cache_SecondCall_ReturnsSameInstance()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var coord = new SweepCoordinator(registry);
        var dim = new SweepDimension.F73Frobenius(
            NValues: new[] { 5 },
            HClasses: new[] { HamiltonianClass.Main },
            ChainOnly: new[] { true });

        var r1 = coord.Sweep(dim);
        var r2 = coord.Sweep(dim);

        Assert.Same(r1, r2);
    }
}
