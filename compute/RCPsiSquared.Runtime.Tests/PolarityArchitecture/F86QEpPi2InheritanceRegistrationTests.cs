using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F86Main;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F86QEpPi2InheritanceRegistrationTests
{
    private const double GammaZero = 0.05;
    private const double GEff = 1.74;

    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            // TPeakLaw's rung-2 edge (2026-06-10) requires AbsorptionTheoremClaim.
            .RegisterAbsorptionTheoremClaim()
            .RegisterF86Main(GammaZero, GEff);

    [Fact]
    public void RegisterF86QEpPi2Inheritance_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF86QEpPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F86QEpPi2Inheritance>());
    }

    [Fact]
    public void RegisterF86QEpPi2Inheritance_AncestorsContainQEpLaw_AndPi2Ladder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF86QEpPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F86QEpPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(QEpLaw), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Fact]
    public void RegisterF86QEpPi2Inheritance_LiveQEpAgreesWithParent()
    {
        // Cross-registry composition check: TwoFactor / GEff = QEpLaw.Value bit-exact.
        var registry = BuildBaseRegistry()
            .RegisterF86QEpPi2Inheritance()
            .Build();

        var inheritance = registry.Get<F86QEpPi2Inheritance>();
        var qEp = registry.Get<QEpLaw>();
        Assert.Equal(qEp.Value, inheritance.LiveQEp, precision: 14);
        Assert.Equal(2.0, inheritance.TwoFactor, precision: 14);
    }

    [Fact]
    public void RegisterF86QEpPi2Inheritance_WithoutPi2Ladder_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                // AbsorptionTheoremClaim registered (TPeakLaw's rung-2 edge, 2026-06-10)
                // so the only missing piece stays the ladder.
                .RegisterAbsorptionTheoremClaim()
                .RegisterF86Main(GammaZero, GEff)
                // Missing: RegisterPi2DyadicLadder
                .RegisterF86QEpPi2Inheritance()
                .Build());
    }
}
