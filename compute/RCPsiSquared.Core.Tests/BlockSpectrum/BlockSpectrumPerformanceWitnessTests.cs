using System.Diagnostics;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Tests.TestHelpers;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

/// <summary>Tests for <see cref="BlockSpectrumPerformanceWitness"/>: claim metadata + a
/// timing assertion at N=5 that the per-block path runs at least
/// <see cref="BlockSpectrumPerformanceWitness.ExpectedSpeedupAtN5"/>× faster than full-L
/// eig with multiset-equal spectrum.
///
/// <para>The timing assertion uses 5 warm-up repeats then 3 measurement repeats per path
/// and asserts on the min of the 3 to be robust against transient JIT/MKL warm-up and CI
/// noise. N=6 timing measurements live in the CLI smoke runs (see
/// <c>compute/RCPsiSquared.Cli/Commands/BlockSpectrumCommand.cs</c>) rather than the unit
/// suite to keep the default test budget cheap.</para></summary>
public class BlockSpectrumPerformanceWitnessTests
{
    // ----------------------------------------------------------------------
    // Claim metadata
    // ----------------------------------------------------------------------

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var claim = MakeClaim();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void DisplayName_IsNonEmpty()
    {
        var claim = MakeClaim();
        Assert.False(string.IsNullOrWhiteSpace(claim.DisplayName));
    }

    [Fact]
    public void Summary_IsNonEmpty()
    {
        var claim = MakeClaim();
        Assert.False(string.IsNullOrWhiteSpace(claim.Summary));
    }

    [Fact]
    public void Constructor_NullParent_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => new BlockSpectrumPerformanceWitness(null!));
    }

    [Fact]
    public void ExpectedSpeedupAtN5_IsAtLeastTwo()
    {
        // The Claim's documented contract: assert >= 2.0× speedup at N=5 to leave headroom
        // for CI machine variance below the observed ~8× margin.
        Assert.True(BlockSpectrumPerformanceWitness.ExpectedSpeedupAtN5 >= 2.0,
            $"ExpectedSpeedupAtN5 = {BlockSpectrumPerformanceWitness.ExpectedSpeedupAtN5} should be >= 2.0 (test contract).");
    }

    // ----------------------------------------------------------------------
    // Timing + correctness witness at N=5
    // ----------------------------------------------------------------------

    [Fact]
    public void BlockEig_AtN5_IsAtLeast2xFasterThanFullL_AndMultisetEqual()
    {
        const int N = 5;
        var L = BuildXYZDephasingL(N, J: 1.0, gamma: 0.5);

        // 5 warm-up repeats: prime JIT, MKL kernels, basis-permutation cache, etc.
        for (int i = 0; i < 5; i++)
        {
            _ = LiouvillianBlockSpectrum.ComputeSpectrum(L, N);
            _ = L.Evd().EigenValues;
        }

        // 3 measurement repeats per path; assert on min to be robust against transient noise.
        // Interleave the two paths so neither benefits from a fresh GC over the other.
        double bestBlockMs = double.PositiveInfinity;
        double bestFullMs = double.PositiveInfinity;
        Complex[] lastBlockEigs = null!;
        Complex[] lastFullEigs = null!;
        var sw = new Stopwatch();
        for (int i = 0; i < 3; i++)
        {
            sw.Restart();
            lastBlockEigs = LiouvillianBlockSpectrum.ComputeSpectrum(L, N);
            sw.Stop();
            double blockMs = sw.Elapsed.TotalMilliseconds;
            if (blockMs < bestBlockMs) bestBlockMs = blockMs;

            sw.Restart();
            lastFullEigs = L.Evd().EigenValues.ToArray();
            sw.Stop();
            double fullMs = sw.Elapsed.TotalMilliseconds;
            if (fullMs < bestFullMs) bestFullMs = fullMs;
        }

        // Correctness: multiset-equal spectrum to 1e-9 tolerance.
        MultisetAssert.NearestNeighbourEqual(lastFullEigs, lastBlockEigs, tolerance: 1e-9, context: $"N={N}");

        // Speedup: block-time * ExpectedSpeedupAtN5 < full-time, i.e. at least
        // ExpectedSpeedupAtN5× faster.
        double speedup = bestFullMs / bestBlockMs;
        Assert.True(
            bestBlockMs * BlockSpectrumPerformanceWitness.ExpectedSpeedupAtN5 < bestFullMs,
            $"Block-eig at N=5 should be >= {BlockSpectrumPerformanceWitness.ExpectedSpeedupAtN5}× faster than full-L. " +
            $"Measured: block min = {bestBlockMs:F2} ms, full min = {bestFullMs:F2} ms, speedup = {speedup:F2}×.");
    }

    // ----------------------------------------------------------------------
    // Helpers
    // ----------------------------------------------------------------------

    private static BlockSpectrumPerformanceWitness MakeClaim()
    {
        var sectors = new JointPopcountSectors();
        var blockSpectrum = new LiouvillianBlockSpectrum(sectors);
        return new BlockSpectrumPerformanceWitness(blockSpectrum);
    }

    private static ComplexMatrix BuildXYZDephasingL(int N, double J, double gamma)
    {
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();
        return PauliDephasingDissipator.BuildZ(H, gammaPerSite);
    }
}
