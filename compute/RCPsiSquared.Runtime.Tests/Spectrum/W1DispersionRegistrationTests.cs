using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Spectrum;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.Spectrum;

namespace RCPsiSquared.Runtime.Tests.Spectrum;

public class W1DispersionRegistrationTests
{
    [Fact]
    public void RegisterW1Dispersion_BuildsOneClaim()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterW1Dispersion(N: 5, J: 1.0, gammaZero: 0.05)
            .Build();

        Assert.Single(registry.All());
        Assert.True(registry.Contains<W1Dispersion>());
    }

    [Fact]
    public void RegisterW1Dispersion_TierIsTier1Derived()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterW1Dispersion(N: 5, J: 1.0, gammaZero: 0.05)
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<W1Dispersion>().Tier);
    }

    [Fact]
    public void RegisterW1Dispersion_AnchorReferencesD10AndAnalyticalSpectrum()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterW1Dispersion(N: 5, J: 1.0, gammaZero: 0.05)
            .Build();

        var anchor = registry.Get<W1Dispersion>().Anchor;
        Assert.Contains("ANALYTICAL_SPECTRUM.md", anchor);
        Assert.Contains("D10_W1_DISPERSION.md", anchor);
    }

    [Fact]
    public void RegisterW1Dispersion_PinnedFrequenciesReproduced()
    {
        // Drift-check: the registered claim returns the same pinned frequency table that
        // experiments/ANALYTICAL_SPECTRUM.md anchors. Same as the Core test, but exercised
        // through the registry to ensure the Schicht-1 wiring does not mutate the values.
        var registry = new ClaimRegistryBuilder()
            .RegisterW1Dispersion(N: 5, J: 1.0, gammaZero: 0.05)
            .Build();

        var d = registry.Get<W1Dispersion>();
        Assert.Equal(0.763932, d.Frequencies[0], precision: 6);
        Assert.Equal(2.763932, d.Frequencies[1], precision: 6);
        Assert.Equal(5.236068, d.Frequencies[2], precision: 6);
        Assert.Equal(7.236068, d.Frequencies[3], precision: 6);
    }

    [Fact]
    public void RegisterW1Dispersion_UniformDecayRateIsTwoGamma()
    {
        // The structural constant the JW track anchors its NEW-NEW Lorentzian width on.
        double γ = 0.05;
        var registry = new ClaimRegistryBuilder()
            .RegisterW1Dispersion(N: 5, J: 1.0, gammaZero: γ)
            .Build();

        Assert.Equal(2.0 * γ, registry.Get<W1Dispersion>().UniformDecayRate, precision: 12);
    }

    [Fact]
    public void RegisterW1Dispersion_NLessThanTwo_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            new ClaimRegistryBuilder()
                .RegisterW1Dispersion(N: 1, J: 1.0, gammaZero: 0.05)
                .Build());
    }
}
