using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F1Pi2InheritanceRegistrationTests
{
    private static ChainSystem DefaultChain() =>
        new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.Heisenberg, Topology: TopologyKind.Chain);

    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterPi2I4MemoryLoop();

    [Fact]
    public void RegisterF1Pi2Inheritance_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF1Pi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F1Pi2Inheritance>());
    }

    [Fact]
    public void RegisterF1Pi2Inheritance_AncestorsContainAllThreeParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF1Pi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F1Pi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(Pi2I4MemoryLoopClaim), ancestors);
    }

    [Fact]
    public void RegisterF1Pi2Inheritance_SignFlipFromZ4_IsMinusOne()
    {
        // Cross-registry verification: F1's "−1" sign flip in "−L" = i² on Z₄
        // memory loop (Layer 1 reading documented in Pi2I4MemoryLoop docstring).
        var registry = BuildBaseRegistry()
            .RegisterF1Pi2Inheritance()
            .Build();

        var f = registry.Get<F1Pi2Inheritance>();
        Assert.Equal(-1.0, f.SignFlipFromZ4.Real, precision: 14);
        Assert.Equal(0.0, f.SignFlipFromZ4.Imaginary, precision: 14);
    }

    [Fact]
    public void RegisterF1Pi2Inheritance_TwoFactorIsTwo()
    {
        var registry = BuildBaseRegistry()
            .RegisterF1Pi2Inheritance()
            .Build();

        Assert.Equal(2.0, registry.Get<F1Pi2Inheritance>().TwoFactor, precision: 14);
    }

    [Fact]
    public void RegisterF1Pi2Inheritance_WithoutF1Family_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterPi2I4MemoryLoop()
                // Missing: RegisterF1Family
                .RegisterF1Pi2Inheritance()
                .Build());
    }

    [Fact]
    public void RegisterF1Pi2Inheritance_WithoutMemoryLoop_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF1Family(DefaultChain())
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                // Missing: RegisterPi2I4MemoryLoop
                .RegisterF1Pi2Inheritance()
                .Build());
    }
}
