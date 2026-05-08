using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.F86Main;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F86LocalGlobalEpLinkRegistrationTests
{
    [Fact]
    public void RegisterF86LocalGlobalEpLink_AddsClaim()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86Main(gammaZero: 0.05, gEff: 1.74)
            .RegisterF86LocalGlobalEpLink()
            .Build();

        Assert.True(registry.Contains<LocalGlobalEpLink>());
    }

    [Fact]
    public void RegisterF86LocalGlobalEpLink_AncestorsContainsChiralAiii()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86Main(gammaZero: 0.05, gEff: 1.74)
            .RegisterF86LocalGlobalEpLink()
            .Build();

        var ancestors = registry.AncestorsOf<LocalGlobalEpLink>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(ChiralAiiiClassification), ancestors);
    }

    [Fact]
    public void RegisterF86LocalGlobalEpLink_TierIsTier2Verified()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86Main(gammaZero: 0.05, gEff: 1.74)
            .RegisterF86LocalGlobalEpLink()
            .Build();

        Assert.Equal(Tier.Tier2Verified, registry.Get<LocalGlobalEpLink>().Tier);
    }

    [Fact]
    public void RegisterF86LocalGlobalEpLink_AnchorReferencesProofF86Qpeak_AndFragileBridge()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86Main(gammaZero: 0.05, gEff: 1.74)
            .RegisterF86LocalGlobalEpLink()
            .Build();

        var anchor = registry.Get<LocalGlobalEpLink>().Anchor;
        Assert.Contains("PROOF_F86_QPEAK.md", anchor);
        Assert.Contains("FRAGILE_BRIDGE.md", anchor);
    }

    [Fact]
    public void RegisterF86LocalGlobalEpLink_WithoutF86Main_Throws_MissingChiralAiiiParent()
    {
        // ChiralAiiiClassification is registered by F86MainRegistration; without it, the
        // edge to ChiralAiiiClassification fails at Build() with MissingParent.
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF86LocalGlobalEpLink()
                .Build());
    }
}
