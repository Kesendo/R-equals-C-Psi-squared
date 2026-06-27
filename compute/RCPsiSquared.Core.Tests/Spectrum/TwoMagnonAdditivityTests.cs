using System.Linq;
using RCPsiSquared.Core.Spectrum;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Spectrum;

/// <summary>The free-fermion-additivity probe for the F89 Door-C CSR sequel
/// (docs/superpowers/plans/2026-06-27-f89-door-c-csr-integrability-sweep.md, review round 1):
/// the sample-size-independent mechanism tracker for "does breaking the (SE,DE) Liouvillian block's
/// free-fermion additivity drive the fixed-q CSR toward Ginibre?". The bra of the (SE,DE) coherence
/// is a two-excitation (two-magnon) state; at Δ=0 the open XX+YY chain is free fermions so the
/// two-magnon energies are EXACT sums of single-magnon energies; the residual measures how fast that
/// additivity breaks as the ZZ-anisotropy Δ turns on (XXZ stays Bethe-integrable, but additivity does
/// NOT survive — that is exactly the distinction the sequel must separate).</summary>
public class TwoMagnonAdditivityTests
{
    private readonly ITestOutputHelper _out;
    public TwoMagnonAdditivityTests(ITestOutputHelper output) => _out = output;

    /// <summary>The load-bearing known-answer anchor: at Δ=0 the open XX+YY chain is free fermions,
    /// so two-magnon energies are EXACT sums of single-magnon energies (no Jordan-Wigner boundary
    /// term on an OPEN chain). The additivity residual must be bit-exact zero.</summary>
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void Residual_AtDeltaZero_IsBitExactZero(int n)
    {
        double r = TwoMagnonAdditivity.Residual(n, j: 1.0, delta: 0.0);
        Assert.True(r < 1e-10, $"N={n}: free-fermion additivity residual {r:E3} should be ~0 at Δ=0");
    }

    /// <summary>Δ turns on the ZZ contact interaction between adjacent magnons, breaking free-fermion
    /// additivity even though the XXZ Hamiltonian stays Bethe-integrable. The residual must be clearly
    /// nonzero — this is the whole point of the probe.</summary>
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void Residual_AtDeltaOne_IsClearlyNonzero(int n)
    {
        double r = TwoMagnonAdditivity.Residual(n, j: 1.0, delta: 1.0);
        Assert.True(r > 1e-3, $"N={n}: additivity residual {r:E3} should be clearly nonzero at Δ=1");
    }

    [Theory]
    [InlineData(4)]
    [InlineData(6)]
    public void SingleMagnonHamiltonian_HasNStates(int n)
    {
        var h = TwoMagnonAdditivity.SingleMagnonHamiltonian(n, j: 1.0, delta: 0.0);
        Assert.Equal(n, h.GetLength(0));
        Assert.Equal(n, h.GetLength(1));
    }

    [Theory]
    [InlineData(4)]  // C(4,2) = 6
    [InlineData(6)]  // C(6,2) = 15
    public void TwoMagnonHamiltonian_HasCN2States(int n)
    {
        int expected = n * (n - 1) / 2;
        var h = TwoMagnonAdditivity.TwoMagnonHamiltonian(n, j: 1.0, delta: 0.0);
        Assert.Equal(expected, h.GetLength(0));
        Assert.Equal(expected, h.GetLength(1));
    }

    /// <summary>Reconnaissance: the additivity residual as Δ turns on, across N. The sample-size-
    /// independent mechanism tracker for the Door-C CSR sequel.</summary>
    [Fact]
    public void Reconnaissance_AdditivityResidualVsDelta()
    {
        _out.WriteLine("  N | Δ=0        | Δ=0.25     | Δ=0.5      | Δ=1.0      | Δ=2.0");
        _out.WriteLine("  --|------------|------------|------------|------------|----------");
        double[] ds = { 0.0, 0.25, 0.5, 1.0, 2.0 };
        foreach (int n in new[] { 5, 6, 7 })
        {
            var vals = ds.Select(d => TwoMagnonAdditivity.Residual(n, 1.0, d)).ToArray();
            _out.WriteLine($"  {n} | " + string.Join(" | ", vals.Select(v => $"{v,10:E3}")));
        }
    }
}
