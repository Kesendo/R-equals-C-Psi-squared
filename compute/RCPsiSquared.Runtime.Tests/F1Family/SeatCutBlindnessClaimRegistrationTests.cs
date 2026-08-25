using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.Tests.F1Family;

/// <summary>Wiring tests for <see cref="SeatCutBlindnessClaim"/> (F157), whose single typed parent
/// is <see cref="F4KernelDimensionByComponentsClaim"/>: F4's extension counts the kernel with the
/// watching on EVERY site and leaves open which seat carries the γ; this claim is the one-seat,
/// popcount-1 corner and closes that question there.</summary>
public class SeatCutBlindnessClaimRegistrationTests
{
    private static ChainSystem DefaultChain(int N = 5) =>
        new ChainSystem(N: N, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.XY, Topology: TopologyKind.Chain);

    private static ClaimRegistry Build() =>
        new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterSeatCutBlindnessClaim()
            .Build();

    [Fact]
    public void RegisterSeatCutBlindness_AddsClaim_Tier1Derived()
    {
        var registry = Build();
        Assert.True(registry.Contains<SeatCutBlindnessClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<SeatCutBlindnessClaim>().Tier);
    }

    [Fact]
    public void RegisterSeatCutBlindness_AncestorsContainTheF4KernelClaim()
    {
        var ancestors = Build().AncestorsOf<SeatCutBlindnessClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F4KernelDimensionByComponentsClaim), ancestors);
    }

    [Fact]
    public void TheParentEdge_IsAlsoHeldAsATypedProperty()
    {
        var registry = Build();
        Assert.Same(registry.Get<F4KernelDimensionByComponentsClaim>(),
                    registry.Get<SeatCutBlindnessClaim>().KernelByComponents);
    }

    [Fact]
    public void TheStatement_KeepsTheOneSeatFence_AndDoesNotInheritTheParentsWording()
    {
        // The parent's law is about UNIFORM dephasing on every site; a repair on 2026-08-24 had to
        // strip a wider "any Σγ > 0" wording from that very claim file. This child must not carry
        // it back in: its hypothesis is one seat, and its sector is popcount 1.
        var claim = Build().Get<SeatCutBlindnessClaim>();
        Assert.Contains("ONE seat", claim.Name);
        Assert.Contains("single-excitation", claim.Name);
        Assert.DoesNotContain("Sigma gamma > 0", claim.Name);
        Assert.DoesNotContain("every site", claim.Name);
    }
}
