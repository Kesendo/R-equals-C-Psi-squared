using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class C2BareDoubledPtfClosedFormRegistrationTests
{
    [Fact]
    public void RegisterC2BareDoubledPtfClosedForm_AddsClaim()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterC2BareDoubledPtfClosedForm()
            .Build();

        Assert.True(registry.Contains<C2BareDoubledPtfClosedForm>());
    }

    [Fact]
    public void RegisterC2BareDoubledPtfClosedForm_TierIsTier1Derived()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterC2BareDoubledPtfClosedForm()
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<C2BareDoubledPtfClosedForm>().Tier);
    }
}
