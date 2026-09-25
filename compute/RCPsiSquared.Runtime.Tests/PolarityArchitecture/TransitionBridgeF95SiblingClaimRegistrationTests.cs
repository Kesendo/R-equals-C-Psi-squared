using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class TransitionBridgeF95SiblingClaimRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterF95AngleAtQuadraticZeroPi2Inheritance();

    [Fact]
    public void RegisteredClaim_ExposesZeroAtBothDiscriminantBoundaries()
    {
        var registry = BuildBaseRegistry()
            .RegisterTransitionBridgeF95SiblingClaim()
            .Build();

        var claim = registry.Get<TransitionBridgeF95SiblingClaim>();
        const double gEff = 0.8;
        double qEp = TransitionBridgeF95SiblingClaim.QEp(gEff);
        Assert.Equal(0.0, claim.CuspAngle(0.25));
        Assert.Equal(0.0, claim.EpClockAngle(1.0, qEp, gEff));
        Assert.Equal(0.0, claim.EpF95Angle(1.0, qEp, gEff));
        Assert.True(claim.EpClockAngleEqualsF95Angle(1.0, qEp, gEff));
    }

    [Fact]
    public void RegisteredClaim_CarriesTheFragileBoundaryWithItsScope()
    {
        var registry = BuildBaseRegistry()
            .RegisterTransitionBridgeF95SiblingClaim()
            .Build();

        var claim = registry.Get<TransitionBridgeF95SiblingClaim>();
        string surface = string.Join("\n", claim.Children.Select(c => $"{c.DisplayName}\n{c.Summary}"));
        Assert.Contains("FRAGILE_BRIDGE boundary", surface, StringComparison.Ordinal);
        Assert.Contains("two qubits per chain", surface, StringComparison.Ordinal);
    }

    [Fact]
    public void RegisteredClaim_HasOnlyTheStandaloneF95Ancestor()
    {
        var registry = BuildBaseRegistry()
            .RegisterTransitionBridgeF95SiblingClaim()
            .Build();

        var ancestors = registry.AncestorsOf<TransitionBridgeF95SiblingClaim>();
        Assert.Single(ancestors);
        Assert.IsType<F95AngleAtQuadraticZeroPi2Inheritance>(ancestors[0]);
        Assert.False(registry.Contains<Pi2DyadicLadderClaim>());
    }
}
