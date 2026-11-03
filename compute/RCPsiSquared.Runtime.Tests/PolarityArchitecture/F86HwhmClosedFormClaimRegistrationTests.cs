using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F86HwhmClosedFormClaimRegistrationTests
{
    [Fact]
    public void RegisterF86HwhmClosedFormClaim_AddsClaim()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86HwhmClosedFormClaim()
            .Build();

        Assert.True(registry.Contains<F86HwhmClosedFormClaim>());
    }

    [Fact]
    public void RegisterF86HwhmClosedFormClaim_TierIsTier2Empirical()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86HwhmClosedFormClaim()
            .Build();

        Assert.Equal(Tier.Tier2Empirical, registry.Get<F86HwhmClosedFormClaim>().Tier);
    }
}
