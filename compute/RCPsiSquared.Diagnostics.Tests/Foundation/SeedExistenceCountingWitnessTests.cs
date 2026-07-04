using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Gate for <see cref="SeedExistenceCountingWitness"/>, the live lab of the
/// seed-existence counting theorem (experiments/F89_SEED_EXISTENCE_REDUCTION.md): on the (1,2)
/// block pencil L(q) = A + qC of the XY chain under uniform Z-dephasing, the real-count identity
/// r(0⁺) − r(∞) = N − 1 for odd N, through the three counting lemmas (N2) n₂ = N − 1,
/// (FF) nullity(C) = ρ, and (N1′) n₆ = 3·Z₃ = ρ (the ordering-sector theorem).
///
/// <para>Two-sided: the zero gates (cross-sector elements, gauge residual) sit beside nonzero
/// controls (the −2/−6 coupling block and, for N ≥ 5, K₆ itself), so a construction or sign
/// error cannot pass silently.</para></summary>
public class SeedExistenceCountingWitnessTests
{
    [Fact]
    public void N5_PinsAllCounts()
    {
        var w = new SeedExistenceCountingWitness(5);
        Assert.Equal(50, w.Dim);
        Assert.Equal(20, w.DimMinus2);
        Assert.Equal(30, w.DimMinus6);
        Assert.Equal(4, w.NullityMinus2);          // (N2) n2 = N - 1
        Assert.Equal(6, w.NullityMinus6);          // (N1') n6 = 3*Z3 = rho
        Assert.Equal(6, w.NullityC);               // (FF)
        Assert.Equal(2, w.ZeroSumTripleCount);     // Z3
        Assert.Equal(6, w.ResonanceCount);         // rho
        Assert.Equal(4, w.Surplus);                // r(0+) - r(inf) = N - 1
        Assert.Equal(0.0, w.CrossSectorMax);       // no-passing, exact
        Assert.Equal(0.0, w.GaugeResidualMax);     // U K_sec U = -H3, exact
        Assert.True(w.CouplingBlockMax > 0.0);     // the -2/-6 coupling is really there
        Assert.True(w.KSixMax > 0.0);              // K6 itself nonzero for N >= 5
    }

    [Fact]
    public void N7_PinsAllCounts()
    {
        var w = new SeedExistenceCountingWitness(7);
        Assert.Equal(147, w.Dim);
        Assert.Equal(6, w.NullityMinus2);
        Assert.Equal(9, w.NullityMinus6);
        Assert.Equal(9, w.NullityC);
        Assert.Equal(3, w.ZeroSumTripleCount);
        Assert.Equal(9, w.ResonanceCount);
        Assert.Equal(6, w.Surplus);
        Assert.Equal(0.0, w.CrossSectorMax);
        Assert.Equal(0.0, w.GaugeResidualMax);
    }

    [Fact]
    public void N3_EdgeCase_KSixIsZeroAndFullKernel()
    {
        // At N = 3 every hop leaves the -6 rung, so K6 = 0 and its kernel is the whole rung.
        var w = new SeedExistenceCountingWitness(3);
        Assert.Equal(3, w.DimMinus6);
        Assert.Equal(3, w.NullityMinus6);
        Assert.Equal(0.0, w.KSixMax);
        Assert.Equal(2, w.Surplus);
        Assert.True(w.CouplingBlockMax > 0.0);     // the coupling is nonzero even at N = 3
    }

    [Theory]
    [InlineData(4)]     // even N is not the claim's regime
    [InlineData(2)]
    [InlineData(11)]    // above the live-build guard
    public void InvalidN_Throws(int n)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new SeedExistenceCountingWitness(n));
    }

    [Fact]
    public void Children_ExposeTheGateNodes()
    {
        var w = new SeedExistenceCountingWitness(5);
        var names = w.Children.Select(c => c.DisplayName).ToList();
        Assert.Contains(names, n => n.Contains("surplus"));
        Assert.Contains(names, n => n.Contains("ordering sectors"));
        Assert.Contains(names, n => n.Contains("open ink"));
    }
}
