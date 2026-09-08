using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class TransitionBridgeF95SiblingClaimRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
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
    public void RegisteredClaim_LeavesFragileThresholdCharacterOpen()
    {
        var registry = BuildBaseRegistry()
            .RegisterTransitionBridgeF95SiblingClaim()
            .Build();

        var claim = registry.Get<TransitionBridgeF95SiblingClaim>();
        string surface = string.Join("\n", claim.Children.Select(c => $"{c.DisplayName}\n{c.Summary}"));
        Assert.Contains("EP/Hopf/Jordan character OPEN", surface, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("FRAGILE_BRIDGE (the EP", surface, StringComparison.OrdinalIgnoreCase);
    }
}
