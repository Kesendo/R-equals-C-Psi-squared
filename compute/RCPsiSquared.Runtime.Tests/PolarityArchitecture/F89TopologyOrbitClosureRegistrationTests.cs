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
    public void RegisterF89_AncestorsContainBothCitedParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89TopologyOrbitClosure()
            .Build();

        var ancestors = registry.AncestorsOf<F89TopologyOrbitClosure>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F73SpatialSumPurityClosurePi2Inheritance), ancestors);
        Assert.Contains(typeof(F71MirrorSymmetryPi2Inheritance), ancestors);
    }
}
