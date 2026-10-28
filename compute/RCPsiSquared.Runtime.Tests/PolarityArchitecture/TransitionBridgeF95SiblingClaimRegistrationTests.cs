using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class TransitionBridgeF95SiblingClaimRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterTwoByTwoExceptionalPointClaim()
            .RegisterF95AngleAtQuadraticZero();

    [Fact]
    public void RegisteredClaim_ExposesZeroAtBothDiscriminantBoundaries()
    {
        var registry = BuildBaseRegistry()
            .RegisterTransitionBridgeF95SiblingClaim()
            .Build();

        var claim = registry.Get<TransitionBridgeF95SiblingClaim>();
        // g_eff = ½ puts Q_EP = 4 on a double, so both routes compute 0.0 there themselves.
        const double gEff = 0.5;
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
    public void RegisteredClaim_HasTheF95LensAndTheOneClosedTwoByTwoAsAncestors()
    {
        var registry = BuildBaseRegistry()
            .RegisterTransitionBridgeF95SiblingClaim()
            .Build();

        var ancestors = registry.AncestorsOf<TransitionBridgeF95SiblingClaim>();
        Assert.Equal(2, ancestors.Count);
        Assert.Contains(ancestors, a => a is F95AngleAtQuadraticZero);
        Assert.Contains(ancestors, a => a is TwoByTwoExceptionalPointClaim);
        Assert.False(registry.Contains<Pi2DyadicLadderClaim>());
    }
}
