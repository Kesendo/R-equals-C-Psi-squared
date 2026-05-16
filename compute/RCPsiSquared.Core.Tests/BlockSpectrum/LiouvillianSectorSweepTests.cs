using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Tests.TestHelpers;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

/// <summary>Witnesses for <see cref="LiouvillianSectorSweep"/>: per-(p_c, p_r) sector
/// dense Evd over the chain XY + Z-dephasing Liouvillian, capped at a user-supplied
/// sectorDim ceiling. Returns the union of computed eigenvalues plus an F1-palindrome
/// witness verifying the spectrum is closed under λ → −λ − 2Σγ within numerical tolerance.
/// At small N where the cap covers every sector the sweep reproduces the full 4^N spectrum
/// bit-for-bit; at larger N some sectors are skipped and the witness operates on the
/// partial-but-symmetric collected set.
///
/// <para>The F1 witness is expected to land at machine precision in this setup (truly XY +
/// Z-dephasing → no F1-Brecher). The break mechanisms (T1 amplitude damping, depolarising
/// noise, transverse-field Hamiltonians — see <see cref="F1.F1OpenQuestions"/> and
/// <c>PalindromeResidualTests.F1_Palindrome_BreaksFor_T1Dissipator</c>) would lift the
/// residual to <c>O(γ)</c> instead of <c>O(N · γ · ε_FP)</c>; this primitive's tight
/// residuals at N=10 confirm no Brecher is silently introduced.</para></summary>
public class LiouvillianSectorSweepTests
{
    private readonly ITestOutputHelper _out;
    public LiouvillianSectorSweepTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void Build_FullCoverage_MatchesFullLSpectrumAsMultiset(int N)
    {
        // At small N with cap = 4^N (= sectorDim upper bound), every sector is included
        // and the collected eigenvalues equal the direct full-L spectrum as a multiset.
        var gamma = Enumerable.Repeat(0.1, N).ToArray();
        int fullCap = 1 << (2 * N);  // 4^N — larger than any sector

        var sweep = LiouvillianSectorSweep.Build(N, gamma, fullCap);

        Assert.Empty(sweep.Skipped);
        Assert.Equal((long)fullCap, sweep.CollectedEigenvalues.LongLength);

        var H = PauliHamiltonian.XYChain(N, 1.0).ToMatrix();
        var fullL = PauliDephasingDissipator.BuildZ(H, gamma);
        var fullEigs = fullL.Evd().EigenValues.ToArray();

        // XY-chain spectra are highly degenerate — sort+zip mis-pairs conjugate clusters.
        // Use greedy nearest-neighbour multiset matching from TestHelpers.
        MultisetAssert.NearestNeighbourEqual(
            sweep.CollectedEigenvalues, fullEigs, tolerance: 1e-9, context: $"N={N}");
        _out.WriteLine($"N={N}: sweep multiset matches full-L spectrum within 1e-9");
    }

    [Theory]
    [InlineData(3, 0.1)]
    [InlineData(4, 0.1)]
    [InlineData(5, 0.05)]
    public void Build_UniformGamma_F1PalindromeHolds(int N, double gammaUniform)
    {
        // F1: for every λ the mirror −λ − 2Σγ is also an eigenvalue. At full coverage
        // under uniform γ the palindrome residual should be at FP-noise level.
        var gamma = Enumerable.Repeat(gammaUniform, N).ToArray();
        int fullCap = 1 << (2 * N);

        var sweep = LiouvillianSectorSweep.Build(N, gamma, fullCap);

        _out.WriteLine($"N={N}, γ={gammaUniform}, Σγ={N * gammaUniform}: " +
                       $"F1 palindrome max residual = {sweep.F1PalindromeResidualMax:G3}");
        Assert.True(sweep.F1PalindromeResidualMax < 1e-9,
            $"F1 palindrome violated: max residual {sweep.F1PalindromeResidualMax:G3}");
    }

    [Fact]
    public void Build_SectorCap_SkipsLargeSectors()
    {
        // At N=5 with cap = 50, the (2,2), (2,3), (3,2), (3,3) sectors should be skipped
        // (their dims are 100, 100, 100, 100). The (1,*) and (0,*) families fit (max dim 25).
        var gamma = Enumerable.Repeat(0.1, 5).ToArray();
        var sweep = LiouvillianSectorSweep.Build(N: 5, gamma, sectorDimCap: 50);

        Assert.NotEmpty(sweep.Skipped);
        foreach (var s in sweep.Skipped)
            Assert.True(s.Dim > 50, $"Skipped sector dim {s.Dim} should exceed cap 50");
        foreach (var s in sweep.Included)
            Assert.True(s.Dim <= 50, $"Included sector dim {s.Dim} should fit cap 50");
        Assert.True(sweep.CoverageFraction < 1.0);
        Assert.True(sweep.CoverageFraction > 0.0);

        _out.WriteLine($"N=5 cap=50: included {sweep.Included.Count}, skipped {sweep.Skipped.Count}, " +
                       $"coverage = {sweep.CoverageFraction:P2}");
    }

    [Fact]
    public void Build_Tier_IsTier1Derived()
    {
        var gamma = Enumerable.Repeat(0.1, 3).ToArray();
        var sweep = LiouvillianSectorSweep.Build(N: 3, gamma, sectorDimCap: 64);
        Assert.Equal(Tier.Tier1Derived, sweep.Tier);
    }

    [Fact]
    public void Build_RejectsGammaWrongLength()
    {
        var wrongGamma = new double[3];
        Assert.Throws<ArgumentException>(() =>
            LiouvillianSectorSweep.Build(N: 5, wrongGamma, sectorDimCap: 100));
    }

    [Theory]
    [InlineData(2500)]   // ~4 % coverage, ~1 min — fast CI smoke
    [InlineData(5500)]   // ~? % coverage, ~few min — deeper reconnaissance, includes (2,3)/(3,2)
                         // and X⊗N partners (dim 5400 each, the next tier of sectors).
    public void Build_AtN10_PartialCoverage_F1PalindromeHolds(int sectorDimCap)
    {
        // F1 palindrome should hold on the symmetric partial collected set at any cap
        // (X⊗N pairing keeps included-sector set closed under the F1 involution).
        var gamma = Enumerable.Repeat(0.05, 10).ToArray();
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var sweep = LiouvillianSectorSweep.Build(N: 10, gamma, sectorDimCap);
        sw.Stop();

        Assert.True(sweep.CollectedEigenvalues.LongLength > 0);
        Assert.True(sweep.F1PalindromeResidualMax < 1e-7,
            $"N=10 partial-spectrum F1 palindrome residual {sweep.F1PalindromeResidualMax:G3} too large");

        _out.WriteLine($"N=10 cap={sectorDimCap}: " +
                       $"included {sweep.Included.Count}, skipped {sweep.Skipped.Count}, " +
                       $"coverage = {sweep.CollectedEigenvalues.LongLength:N0} / {sweep.FullDim:N0} " +
                       $"= {sweep.CoverageFraction:P3}, wall = {sw.Elapsed.TotalSeconds:F1} s");
        _out.WriteLine($"  F1 palindrome max residual = {sweep.F1PalindromeResidualMax:G3}");
    }
}
