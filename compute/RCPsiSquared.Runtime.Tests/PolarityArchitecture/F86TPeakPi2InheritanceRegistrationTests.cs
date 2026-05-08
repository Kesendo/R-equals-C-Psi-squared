using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F86Main;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F86TPeakPi2InheritanceRegistrationTests
{
    private const double GammaZero = 0.05;
    private const double GEff = 1.74;

    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF86Main(GammaZero, GEff);

    [Fact]
    public void RegisterF86TPeakPi2Inheritance_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF86TPeakPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F86TPeakPi2Inheritance>());
    }

    [Fact]
    public void RegisterF86TPeakPi2Inheritance_AncestorsContainTPeakLaw_AndPi2Ladder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF86TPeakPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F86TPeakPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(TPeakLaw), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Fact]
    public void RegisterF86TPeakPi2Inheritance_LiveTPeakAgreesWithParent()
    {
        // Cross-registry composition check: 1/(FourFactor · GammaZero) = TPeakLaw.Value
        // bit-exact. The "4" denominator is a_{-1} = d² for N=1.
        var registry = BuildBaseRegistry()
            .RegisterF86TPeakPi2Inheritance()
            .Build();

        var inheritance = registry.Get<F86TPeakPi2Inheritance>();
        var tPeak = registry.Get<TPeakLaw>();
        Assert.Equal(tPeak.Value, inheritance.LiveTPeak, precision: 14);
        Assert.Equal(4.0, inheritance.FourFactor, precision: 14);
        Assert.Equal(0.25, inheritance.OneOverFourFactor, precision: 14);
    }

    [Fact]
    public void RegisterF86TPeakPi2Inheritance_MirrorPartnerProductIsOne()
    {
        // The inversion symmetry a_{-1} · a_3 = 4 · 1/4 = 1; both readings of the same fact.
        var registry = BuildBaseRegistry()
            .RegisterF86TPeakPi2Inheritance()
            .Build();

        var inheritance = registry.Get<F86TPeakPi2Inheritance>();
        Assert.Equal(1.0, inheritance.FourFactor * inheritance.OneOverFourFactor, precision: 14);
    }
}
