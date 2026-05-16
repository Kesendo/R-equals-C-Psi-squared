using System.Numerics;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Tests.TestHelpers;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

/// <summary>Cross-validation of the BlockSpectrum primitives against the EQ-014 ground
/// truth — the full N=7 Liouvillian eigendecomposition produced by
/// <c>RCPsiSquared.Compute.Program.RunPtfExport</c> (~48 h wall time, ~16 GB peak memory).
///
/// <para>Existing BlockSpectrum tests verify per-block eig vs direct full-L only up to
/// N=5 (direct full-L at N=7 needs 8 GB just to materialise the matrix). The EQ-014
/// precomputed data lets us extend the cross-validation to N=7 without re-running the
/// full-L diagonalisation each test run.</para>
///
/// <para>Tests skip gracefully when the precomputed file is absent — only the developer
/// who ran the PTF export has it locally, and CI environments without that data should
/// not fail.</para>
/// </summary>
public class Eq014N7CrossValidationTests
{
    private readonly ITestOutputHelper _out;
    public Eq014N7CrossValidationTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Eq014_EigenvalueFile_ExistsAndHasExpectedShape()
    {
        if (!Eq014GroundTruth.IsAvailable(Eq014GroundTruth.EigenvaluesFileName))
        {
            _out.WriteLine($"Skipping: {Eq014GroundTruth.EigenvaluesFileName} not present " +
                           $"at {Eq014GroundTruth.ResultsDirectory}. " +
                           $"Generate via `cd compute/RCPsiSquared.Compute && dotnet run -c Release -- ptf`.");
            return;
        }

        var eigs = Eq014GroundTruth.LoadEigenvalues();
        Assert.Equal(Eq014GroundTruth.LiouvilleDim, eigs.Length);

        // Structural sanity from the metadata:
        // - stationary count: F4 says N+1 = 8 stationary modes (|λ| < 1e-10)
        // - palindrome center: -Σγ = -7·0.05 = -0.35
        int stationary = eigs.Count(z => z.Magnitude < 1e-10);
        Assert.Equal(8, stationary);

        double sumGamma = Eq014GroundTruth.N * Eq014GroundTruth.Gamma;
        double palindromeCenter = -sumGamma;
        // Spectrum mean of real parts = palindrome center (F1: λ pairs with −λ−2Σγ around −Σγ).
        double meanRe = eigs.Sum(z => z.Real) / eigs.Length;
        Assert.Equal(palindromeCenter, meanRe, precision: 6);

        _out.WriteLine($"eq014 N=7: {eigs.Length} eigenvalues, {stationary} stationary, " +
                       $"mean Re(λ) = {meanRe:F6} (palindrome center = {palindromeCenter:F6})");
    }

    [Fact]
    public void LiouvillianBlockSpectrum_ComputeSpectrumPerBlock_AtN7_MatchesEq014()
    {
        if (!Eq014GroundTruth.IsAvailable(Eq014GroundTruth.EigenvaluesFileName))
        {
            _out.WriteLine($"Skipping: {Eq014GroundTruth.EigenvaluesFileName} not present.");
            return;
        }

        var groundTruth = Eq014GroundTruth.LoadEigenvalues();
        var H = PauliHamiltonian.XYChain(Eq014GroundTruth.N, Eq014GroundTruth.J).ToMatrix();
        var gamma = Enumerable.Repeat(Eq014GroundTruth.Gamma, Eq014GroundTruth.N).ToArray();

        var sw = System.Diagnostics.Stopwatch.StartNew();
        var blockSpectrum = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gamma, Eq014GroundTruth.N);
        sw.Stop();

        Assert.Equal(Eq014GroundTruth.LiouvilleDim, blockSpectrum.Length);

        var matchSw = System.Diagnostics.Stopwatch.StartNew();
        MultisetAssert.NearestNeighbourEqual(
            blockSpectrum, groundTruth, tolerance: 1e-9,
            context: $"N=7 ComputeSpectrumPerBlock vs eq014");
        matchSw.Stop();

        _out.WriteLine($"N=7 per-block eig: {sw.Elapsed.TotalSeconds:F1} s; " +
                       $"multiset match: {matchSw.Elapsed.TotalSeconds:F1} s; " +
                       $"all {blockSpectrum.Length} eigenvalues match eq014 within 1e-9.");
    }

    [Fact]
    public void LiouvillianSectorSweep_AtN7_FullCoverage_MatchesEq014()
    {
        if (!Eq014GroundTruth.IsAvailable(Eq014GroundTruth.EigenvaluesFileName))
        {
            _out.WriteLine($"Skipping: {Eq014GroundTruth.EigenvaluesFileName} not present.");
            return;
        }

        var groundTruth = Eq014GroundTruth.LoadEigenvalues();
        var gamma = Enumerable.Repeat(Eq014GroundTruth.Gamma, Eq014GroundTruth.N).ToArray();

        // Max sector dim at N=7 = C(7,3)² = 35² = 1225; cap 1500 covers every sector.
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var sweep = LiouvillianSectorSweep.Build(Eq014GroundTruth.N, gamma, sectorDimCap: 1500);
        sw.Stop();

        Assert.Empty(sweep.Skipped);
        Assert.Equal(Eq014GroundTruth.LiouvilleDim, (int)sweep.CollectedEigenvalues.LongLength);

        var matchSw = System.Diagnostics.Stopwatch.StartNew();
        MultisetAssert.NearestNeighbourEqual(
            sweep.CollectedEigenvalues, groundTruth, tolerance: 1e-9,
            context: $"N=7 LiouvillianSectorSweep vs eq014");
        matchSw.Stop();

        _out.WriteLine($"N=7 sector sweep (cap 1500): {sw.Elapsed.TotalSeconds:F1} s; " +
                       $"multiset match: {matchSw.Elapsed.TotalSeconds:F1} s; " +
                       $"F1 palindrome residual = {sweep.F1PalindromeResidualMax:G3}; " +
                       $"all {sweep.CollectedEigenvalues.LongLength:N0} eigenvalues match eq014 within 1e-9.");
    }
}
