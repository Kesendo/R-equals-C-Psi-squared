using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F89TopologyOrbitClosureRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterHalfIntegerMirror(N: 5)
            .RegisterF70DeltaNSelectionRulePi2Inheritance()
            .RegisterF72BlockDiagonalPurityPi2Inheritance()
            .RegisterF73SpatialSumPurityClosurePi2Inheritance()
            .RegisterF71MirrorSymmetryPi2Inheritance();

    [Fact]
    public void RegisterF89_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89TopologyOrbitClosure()
            .Build();

        Assert.True(registry.Contains<F89TopologyOrbitClosure>());
    }

    [Fact]
    public void RegisterF89_AncestorsContainAllParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89TopologyOrbitClosure()
            .Build();

        var ancestors = registry.AncestorsOf<F89TopologyOrbitClosure>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(F73SpatialSumPurityClosurePi2Inheritance), ancestors);
        Assert.Contains(typeof(F71MirrorSymmetryPi2Inheritance), ancestors);
    }

    [Fact]
    public void RegisterF89_DecayCoefficient_MatchesPi2LadderAnchor()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89TopologyOrbitClosure()
            .Build();

        var f89 = registry.Get<F89TopologyOrbitClosure>();
        Assert.Equal(4.0, f89.DecayRateCoefficient, precision: 14);
        Assert.Equal(4.0, f89.OscillationFrequencyCoefficient, precision: 14);
        Assert.True(f89.Pi2DoublingConsistent());
    }
}
