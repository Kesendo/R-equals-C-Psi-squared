using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F44AlgebraicPairRateLogIdentityPi2InheritanceRegistrationTests
{
    private static ChainSystem DefaultChain() =>
        new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.Heisenberg, Topology: TopologyKind.Chain);

    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterPi2I4MemoryLoop()
            // Welle 7 (2026-05-25) wired F1 → F61 chain; expand deps so
            // F1Pi2Inheritance resolves. See F1Pi2InheritanceRegistrationTests.
            .RegisterF88bPopcountCoherence()
            .RegisterF88bStaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror()
            .RegisterF38BitAInvolutionInheritance()
            .RegisterF63BitAReference()
            .RegisterF38Pi2InvolutionPi2Inheritance()
            .RegisterF63LCommutesPi2Pi2Inheritance()
            .RegisterF61BitAParityPi2Inheritance()
            .RegisterF1Pi2Inheritance();

    [Fact]
    public void RegisterF44_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF44AlgebraicPairRateLogIdentityPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F44AlgebraicPairRateLogIdentityPi2Inheritance>());
    }

    [Fact]
    public void RegisterF44_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF44AlgebraicPairRateLogIdentityPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F44AlgebraicPairRateLogIdentityPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF44_AncestorsContainF1AndPi2DyadicLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF44AlgebraicPairRateLogIdentityPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F44AlgebraicPairRateLogIdentityPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1Pi2Inheritance), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(Pi2I4MemoryLoopClaim), ancestors);
    }

    [Theory]
    // F68 N=3 example: d_fast = 0.075, d_slow = 0.025, Δd = 0.05, Σγ = 0.05
    // → ln(d_fast/d_slow) = ln(3) ≈ 1.0986
    [InlineData(0.05, 0.05, 1.0986122886681098)]
    [InlineData(0.0, 0.05, 0.0)]
    public void RegisterF44_LogRatioAcrossRegistry(double deltaD, double totalGamma, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF44AlgebraicPairRateLogIdentityPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F44AlgebraicPairRateLogIdentityPi2Inheritance>().LogRatio(deltaD, totalGamma), precision: 12);
    }

    [Fact]
    public void RegisterF44_PalindromicSumHoldsForF68Example()
    {
        var registry = BuildBaseRegistry()
            .RegisterF44AlgebraicPairRateLogIdentityPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F44AlgebraicPairRateLogIdentityPi2Inheritance>().PalindromicSumHolds(
            dFast: 0.075, dSlow: 0.025, totalGamma: 0.05));
    }

    [Fact]
    public void RegisterF44_RejectsValuesOutsideFiniteOrderedPairDomain()
    {
        var claim = BuildBaseRegistry()
            .RegisterF44AlgebraicPairRateLogIdentityPi2Inheritance()
            .Build()
            .Get<F44AlgebraicPairRateLogIdentityPi2Inheritance>();

        Assert.Throws<ArgumentOutOfRangeException>(() => claim.LogRatio(-0.01, 0.05));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.LogRatioLinearApproximation(0.10, 0.05));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.LogRatio(double.NaN, 0.05));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.LogRatioLinearApproximation(0.01, double.PositiveInfinity));
    }

    [Fact]
    public void RegisterF44_ExposesOnlyTheAlgebraicInverseRateCoefficient()
    {
        var registry = BuildBaseRegistry()
            .RegisterF44AlgebraicPairRateLogIdentityPi2Inheritance()
            .Build();

        var claim = registry.Get<F44AlgebraicPairRateLogIdentityPi2Inheritance>();
        Assert.Equal(20.0, claim.InverseRateLinearCoefficient(0.05), precision: 12);
        Assert.DoesNotContain(claim.GetType().GetMembers(), member =>
            member.Name.Contains("Temperature", StringComparison.OrdinalIgnoreCase) ||
            member.Name.Contains("Crooks", StringComparison.OrdinalIgnoreCase) ||
            member.Name.Contains("Jarzynski", StringComparison.OrdinalIgnoreCase));
    }

    [Fact]
    public void RegisterF44_WithoutF1_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF1Family(DefaultChain())
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterPi2I4MemoryLoop()
                // Missing: RegisterF1Pi2Inheritance
                .RegisterF44AlgebraicPairRateLogIdentityPi2Inheritance()
                .Build());
    }
}
