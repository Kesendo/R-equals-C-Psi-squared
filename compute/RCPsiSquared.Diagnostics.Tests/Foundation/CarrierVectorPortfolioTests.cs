using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class CarrierVectorPortfolioTests
{
    private readonly Xunit.Abstractions.ITestOutputHelper _out;
    public CarrierVectorPortfolioTests(Xunit.Abstractions.ITestOutputHelper output) => _out = output;

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
        // −Re(λ) = 2·Σ_x γ_x·⟨Δ_x⟩, reproducing sip_carrier_channels.py portfolio rows.
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

    private static Matrix<Complex> Heisenberg2(double J)
    {
        var X = Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 0, 1 }, { 1, 0 } });
        var Y = Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 0, -Complex.ImaginaryOne }, { Complex.ImaginaryOne, 0 } });
        var Z = Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 1, 0 }, { 0, -1 } });
        return (X.KroneckerProduct(X) + Y.KroneckerProduct(Y) + Z.KroneckerProduct(Z)).Multiply((Complex)J);
    }

    // Heisenberg chain J*sum_l (X_lX_{l+1}+Y_lY_{l+1}+Z_lZ_{l+1}), number-conserving (block path applies).
    private static Matrix<Complex> HeisChain(int n, double J)
    {
        var X = Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 0, 1 }, { 1, 0 } });
        var Y = Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 0, -Complex.ImaginaryOne }, { Complex.ImaginaryOne, 0 } });
        var Z = Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 1, 0 }, { 0, -1 } });
        var I2 = Matrix<Complex>.Build.DenseIdentity(2);
        var H = Matrix<Complex>.Build.Dense(1 << n, 1 << n);
        for (int l = 0; l < n - 1; l++)
            foreach (var P in new[] { X, Y, Z })
            {
                Matrix<Complex>? term = null;
                for (int s = 0; s < n; s++)
                {
                    var op = (s == l || s == l + 1) ? P : I2;
                    term = term is null ? op : term.KroneckerProduct(op);
                }
                H += term!.Multiply((Complex)J);
            }
        return H;
    }

    [Fact]
    public void DecomposeBlocked_MatchesFullDecompose_ForNumberConservingH()
    {
        // Non-uniform carrier (the case where the site-labeling convention actually matters).
        var channels = new[] { new ChannelRate("q0", 0.05), new ChannelRate("q1", 0.13), new ChannelRate("q2", 0.20) };
        var H = HeisChain(3, 1.0);
        Assert.True(CarrierVectorPortfolio.IsNumberConserving(H), "Heisenberg chain should be number-conserving");

        var full = CarrierVectorPortfolio.Decompose(3, H, channels);
        var blocked = CarrierVectorPortfolio.DecomposeBlocked(3, H, channels);

        // Same count of modes (4^3 = 64) and the same spectrum (sorted rates, bit-exact).
        Assert.Equal(full.Modes.Count, blocked.Modes.Count);
        Assert.Equal(64, blocked.Modes.Count);
        var fullRates = full.Modes.Select(m => Math.Round(m.ActualDecayRate, 9)).OrderBy(x => x).ToArray();
        var blockedRates = blocked.Modes.Select(m => Math.Round(m.ActualDecayRate, 9)).OrderBy(x => x).ToArray();
        Assert.Equal(fullRates, blockedRates);

        // The blocked portfolios satisfy the carrier-vector law: if the reversed-gamma / (N-1-s)-bit
        // convention were wrong, this residual would be O(1), not machine zero.
        Assert.True(blocked.MaxLawResidual < 1e-9, $"blocked law residual {blocked.MaxLawResidual:E3} exceeds 1e-9");
    }

    [Fact]
    public void IsNumberConserving_SeesTheDifference()
    {
        Assert.True(CarrierVectorPortfolio.IsNumberConserving(HeisChain(3, 1.0)));
        // sum_l X_l Y_{l+1} changes popcount by 0, +-2: not number-conserving.
        var X = Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 0, 1 }, { 1, 0 } });
        var Y = Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 0, -Complex.ImaginaryOne }, { Complex.ImaginaryOne, 0 } });
        var I2 = Matrix<Complex>.Build.DenseIdentity(2);
        var xy = X.KroneckerProduct(Y).KroneckerProduct(I2) + I2.KroneckerProduct(X).KroneckerProduct(Y);
        Assert.False(CarrierVectorPortfolio.IsNumberConserving(xy));
    }

    [Fact]
    public void Decompose_FromLiveLiouvillian_ClosesTheLawAndReadsSubsetSums()
    {
        var channels = new[] { new ChannelRate("a", 0.05), new ChannelRate("b", 0.20) };

        // Coupled Heisenberg: each mode's portfolio, read from the eigenvector, reproduces
        // its own decay rate. The law closes on a live, non-trivial Liouvillian.
        var coupled = CarrierVectorPortfolio.Decompose(2, Heisenberg2(1.0), channels);
        Assert.Equal(16, coupled.Modes.Count);
        Assert.True(coupled.MaxLawResidual < 1e-9,
            $"coupled law residual {coupled.MaxLawResidual:E3} exceeds 1e-9");

        // H = 0: pure dephasing, rates are the exact subset-sums of {2γ_a, 2γ_b} = {0, .1, .4, .5}.
        var zero = Matrix<Complex>.Build.Dense(4, 4);
        var deph = CarrierVectorPortfolio.Decompose(2, zero, channels);
        var rates = deph.Modes.Select(m => Math.Round(m.ActualDecayRate, 6)).Distinct().ToList();
        Assert.Equal(4, rates.Count);
        foreach (var expected in new[] { 0.0, 0.1, 0.4, 0.5 })
            Assert.Contains(rates, r => Math.Abs(r - expected) < 1e-6);
        Assert.True(deph.MaxLawResidual < 1e-9);
    }

    [Fact]
    public void Look_GoodQubitProtectionLeaksWhenCoupled()
    {
        // A slow "good" qubit (γ=0.01) and a fast "bad" one (γ=1.0).
        var channels = new[] { new ChannelRate("good", 0.01), new ChannelRate("bad", 1.0) };
        var coupled = CarrierVectorPortfolio.Decompose(2, Heisenberg2(1.0), channels);
        var decoupled = CarrierVectorPortfolio.Decompose(2, Matrix<Complex>.Build.Dense(4, 4), channels);

        void Show(string label, CarrierPortfolioSpectrum s)
        {
            _out.WriteLine($"{label}  (slowest nonzero rate = {s.SlowestRate:F4})");
            foreach (var m in s.Modes.Where(m => m.ActualDecayRate > 1e-9)
                                     .OrderBy(m => m.ActualDecayRate).Take(2))
            {
                var p = string.Join("  ", m.Portfolio.Activity.Select(a => $"{a.Channel} {100 * a.Delta:F0}%"));
                _out.WriteLine($"    rate {m.ActualDecayRate:F4}   {p}");
            }
        }
        Show("DECOUPLED (good qubit isolated)", decoupled);
        Show("COUPLED   (good qubit + bad neighbour)", coupled);

        // Decoupled: the good qubit's coherence is protected at exactly 2·γ_good.
        Assert.Equal(0.02, decoupled.SlowestRate, precision: 6);
        // Coupled: the coupling leaks the bad channel in and raises the protected floor.
        Assert.True(coupled.SlowestRate > decoupled.SlowestRate + 1e-6,
            $"expected leak: coupled {coupled.SlowestRate:F4} should exceed decoupled {decoupled.SlowestRate:F4}");
    }
}
