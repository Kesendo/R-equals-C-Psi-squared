using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.F86Main;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.Tests.F86Main;

public class F86MainRegistrationTests
{
    [Fact]
    public void RegisterF86Main_BuildsFiveClaims()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86Main(gammaZero: 0.05, gEff: 1.0)
            .Build();

        Assert.Equal(5, registry.All().Count());
        Assert.True(registry.Contains<ChiralAiiiClassification>());
        Assert.True(registry.Contains<DressedModeWeightClaim>());
        Assert.True(registry.Contains<F71MirrorInvariance>());
        Assert.True(registry.Contains<TPeakLaw>());
        Assert.True(registry.Contains<QEpLaw>());
    }

    [Fact]
    public void RegisterF86Main_AllTier1Derived()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86Main(gammaZero: 0.05, gEff: 1.0)
            .Build();

        Assert.All(registry.All(), c => Assert.Equal(Tier.Tier1Derived, c.Tier));
    }

    [Fact]
    public void RegisterF86Main_TPeakLaw_UsesProvidedGamma()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86Main(gammaZero: 0.1, gEff: 1.0)
            .Build();

        var tpeak = registry.Get<TPeakLaw>();
        // t_peak = 1/(4γ₀) = 1/0.4 = 2.5
        Assert.Equal(2.5, tpeak.Value, precision: 10);
    }

    [Fact]
    public void RegisterF86Main_QEpLaw_UsesProvidedGEff()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86Main(gammaZero: 0.05, gEff: 4.0)
            .Build();

        var qep = registry.Get<QEpLaw>();
        // Q_EP = 2/g_eff = 2/4 = 0.5
        Assert.Equal(0.5, qep.Value, precision: 10);
    }
}
