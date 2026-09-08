using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.F86Main;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F86LocalGlobalEpLinkRegistrationTests
{
    // TPeakLaw's rung-2 edge (2026-06-10) requires the absorption chain under
    // RegisterF86Main: Pi2Family + Pi2DyadicLadder + AbsorptionTheoremClaim.
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterF1PalindromeIdentity()
            .RegisterF86Main(gammaZero: 0.05, gEff: 1.74);

    [Fact]
    public void RegisterF86LocalGlobalEpLink_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF86LocalGlobalEpLink()
            .Build();

        Assert.True(registry.Contains<LocalGlobalEpLink>());
    }

    [Fact]
    public void RegisterF86LocalGlobalEpLink_AncestorsContainsShiftedGeneratorPTypeClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF86LocalGlobalEpLink()
            .Build();

        var ancestors = registry.AncestorsOf<LocalGlobalEpLink>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(ShiftedGeneratorSectorwisePClaim), ancestors);
    }

    [Fact]
    public void RegisterF86LocalGlobalEpLink_TierIsOpenQuestion()
    {
        // Demoted Tier2Verified → OpenQuestion (F86a repair): the coarse Q grid did not sample
        // the narrow real-axis coalescences later certified at N=5,7,9, and the prior magnitudes
        // are grid-sensitive. The Tier1Derived shifted-generator P-type parent
        // (strength 5) still dominates OpenQuestion (strength 1), so the edge survives.
        var registry = BuildBaseRegistry()
            .RegisterF86LocalGlobalEpLink()
            .Build();

        Assert.Equal(Tier.OpenQuestion, registry.Get<LocalGlobalEpLink>().Tier);
    }

    [Fact]
    public void RegisterF86LocalGlobalEpLink_AnchorReferencesCurrentEpMechanism_AndFragileBridge()
    {
        var registry = BuildBaseRegistry()
            .RegisterF86LocalGlobalEpLink()
            .Build();

        var anchor = registry.Get<LocalGlobalEpLink>().Anchor;
        Assert.Contains("PROOF_F86A_EP_MECHANISM.md", anchor);
        Assert.Contains("FRAGILE_BRIDGE.md", anchor);
    }

    [Fact]
    public void RegisterF86LocalGlobalEpLink_WithoutF86Main_Throws_MissingPTypeParent()
    {
        // ShiftedGeneratorSectorwisePClaim is registered by F86MainRegistration; without it,
        // the parent edge fails at Build() with MissingParent.
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF86LocalGlobalEpLink()
                .Build());
    }
}
