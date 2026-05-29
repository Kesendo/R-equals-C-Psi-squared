using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class CarrierVectorPortfolioTests
{
    // The Si:P channel stack: clocks spanning decades (charge fast, nuclear protected).
    private static CarrierVector SiPCarrier() => new(new[]
    {
        new ChannelRate("e", 1.0),    // electron spin, fast
        new ChannelRate("n", 0.001),  // nuclear spin, very slow / protected
        new ChannelRate("c", 10.0),   // charge, sensitive
        new ChannelRate("v", 0.1),    // valley, medium
    });

    private static ChannelDifferencePortfolio Portfolio(params (string Channel, double Delta)[] a) =>
        new(a.Select(x => new ChannelActivity(x.Channel, x.Delta)).ToList());

    [Fact]
    public void DecayRate_SiPPortfolios_MatchWeightedSum()
    {
        // −Re(λ) = 2·Σ_x γ_x·⟨Δ_x⟩, reproducing _sip_carrier_channels.py portfolio rows.
        var g = SiPCarrier();
        Assert.Equal(2.002, Portfolio(("e", 1), ("n", 1), ("c", 0), ("v", 0)).DecayRate(g), precision: 9);
        Assert.Equal(0.200, Portfolio(("e", 0), ("n", 0), ("c", 0), ("v", 1)).DecayRate(g), precision: 9);
        Assert.Equal(21.201, Portfolio(("e", 0.5), ("n", 0.5), ("c", 1), ("v", 1)).DecayRate(g), precision: 9);
    }

    [Fact]
    public void NuclearOnlyMode_IsProtected_AndSitsAtTheSlowFloor()
    {
        var g = SiPCarrier();
        var pureNuclear = Portfolio(("e", 0), ("n", 1), ("c", 0), ("v", 0));
        // Sweet-spot condition: zero difference stored in the fast channels.
        Assert.True(pureNuclear.IsProtectedFrom(new[] { "e", "c", "v" }));
        // Then only the slow nuclear clock remains: 2·γ_n.
        Assert.Equal(0.002, pureNuclear.DecayRate(g), precision: 9);
    }

    [Fact]
    public void DominantChannel_IsTheHeaviestWeightedActivity()
    {
        // A 50/50 electron-nuclear mode loses its coherence through the electron,
        // because γ_e ≫ γ_n even at equal stored difference.
        var g = SiPCarrier();
        Assert.Equal("e", Portfolio(("e", 0.5), ("n", 0.5)).DominantChannel(g));
    }

    [Fact]
    public void UniformCarrier_ReducesToPopcount()
    {
        var g = new CarrierVector(new[]
        {
            new ChannelRate("0", 0.1), new ChannelRate("1", 0.1), new ChannelRate("2", 0.1),
        });
        Assert.True(g.IsUniform());
        // A sharp two-channel-difference mode decays at 2γ·popcount = 2·0.1·2.
        var depthTwo = Portfolio(("0", 1), ("1", 1), ("2", 0));
        Assert.Equal(0.4, depthTwo.DecayRate(g), precision: 9);
    }

    [Fact]
    public void Spectrum_HoldsAllModesAtOnce_WithLawResidualZero()
    {
        var g = SiPCarrier();
        var modes = new List<CarrierMode>
        {
            new(0.002, Portfolio(("e", 0), ("n", 1), ("c", 0), ("v", 0))),
            new(2.002, Portfolio(("e", 1), ("n", 1), ("c", 0), ("v", 0))),
            new(0.200, Portfolio(("e", 0), ("n", 0), ("c", 0), ("v", 1))),
        };
        var spec = new CarrierPortfolioSpectrum(g, modes);
        Assert.Equal(0.002, spec.SlowestRate, precision: 9);
        Assert.Equal(2.002, spec.FastestRate, precision: 9);
        Assert.Equal(0.002 + 2.002 + 0.200, spec.Budget, precision: 9);
        // Portfolios are the true activities, so the law holds bit-exact across the set.
        Assert.True(spec.MaxLawResidual < 1e-9,
            $"law residual {spec.MaxLawResidual:E3} exceeds 1e-9");
    }
}
