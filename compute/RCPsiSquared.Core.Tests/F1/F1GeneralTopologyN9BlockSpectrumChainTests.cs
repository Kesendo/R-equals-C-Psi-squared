using System.Diagnostics;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Tests.TestHelpers;
using Xunit;
using Xunit.Abstractions;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.F1;

/// <summary>F1 general-topology N=9 chain spot-check via
/// <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/>.
///
/// <para><b>LP64 bridge landed 2026-05-18.</b> At N=9 the max joint-popcount block is
/// C(9, 4) · C(9, 5) = 15 876, i.e. 15 876² × 16 B ≈ 4 GB, exceeding MathNet's LP64
/// <c>MklLinearAlgebraProvider.EigenDecomp</c> 2 GB single-native-array marshalling ceiling.
/// <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/> now auto-routes blocks
/// larger than <see cref="LiouvillianBlockSpectrum.Lp64ComplexCeiling"/> (11 585²) through
/// <c>RCPsiSquared.Core.Numerics.MklDirect</c>'s NativeMemory + ILP64-aware LAPACK path,
/// bypassing the marshaller cap. Bit-exact parity vs MathNet at small N is witnessed by
/// <c>PerBlockLiouvillianBuilderNativeMemoryParityTests</c>.</para>
///
/// <para>The test is opt-in under <c>[Trait("Category", "SLOW_N9")]</c> because the largest
/// sector pair (C(9, 4) · C(9, 5) = 15 876² each) holds ~4 GB native memory per block during
/// its zgeev call; the per-block serialisation built into the MklDirect branch
/// (<c>EigenPath.MklDirectNative</c>) caps wall-time but stays well inside the dev machine's
/// 128 GB envelope.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md</c> (scale-frontier section
/// updated with the bridge landing), <see cref="F1GeneralTopologyVerifiedClaim"/>,
/// <see cref="F1SpectrumStatistics"/> (shared metrics utility),
/// <see cref="LiouvillianBlockSpectrum"/> (the per-block dispatch + Lp64ComplexCeiling
/// constant), <c>compute/RCPsiSquared.Core/Numerics/MklDirect.cs</c> (the NativeMemory +
/// ILP64 LAPACK route that bridges past LP64 marshalling).</para></summary>
public class F1GeneralTopologyN9BlockSpectrumChainTests
{
    private readonly ITestOutputHelper _out;

    public F1GeneralTopologyN9BlockSpectrumChainTests(ITestOutputHelper output) => _out = output;

    /// <summary>Hard ceiling on per-block Evd dimension for MathNet's LP64 MKL provider.
    /// The native marshaller <c>MngdNativeArrayMarshaler.ConvertSpaceToNative</c> caps
    /// individual <c>Complex[]</c> P/Invoke arrays at 2 GB; with 16 B per Complex that is
    /// 134 217 728 elements, i.e. a 11 585² square matrix. Block size 15 876² (N=9 max)
    /// exceeds this; block size 4 900² (N=8 max) is comfortably below.
    /// <para>Kept in lockstep with <see cref="LiouvillianBlockSpectrum.Lp64ComplexCeiling"/>
    /// (the central constant the production dispatch uses) so a downstream reader sees the
    /// same number from either side.</para></summary>
    public const int Lp64EvdSquareMatrixCeiling = LiouvillianBlockSpectrum.Lp64ComplexCeiling;

    private static Bond[] ChainBonds(int N) =>
        Enumerable.Range(0, N - 1).Select(i => new Bond(i, i + 1, 1.0)).ToArray();

    /// <summary>Build H = (J/4) Σ_b (X_i X_j + Y_i Y_j + Z_i Z_j) on the chain.
    /// Popcount-conserving so it stays inside the
    /// <see cref="JointPopcountSectorBuilder"/> block-infrastructure domain.</summary>
    private static ComplexMatrix BuildHeisenbergGraphHamiltonian(int N, IReadOnlyList<Bond> bonds, double J = 1.0)
    {
        var terms = new (PauliLetter, PauliLetter, Complex)[]
        {
            (PauliLetter.X, PauliLetter.X, J / 4.0),
            (PauliLetter.Y, PauliLetter.Y, J / 4.0),
            (PauliLetter.Z, PauliLetter.Z, J / 4.0),
        };
        return PauliHamiltonian.Bilinear(N, bonds, terms).ToMatrix();
    }

    [Fact]
    [Trait("Category", "SLOW_N9")]
    public void Chain_HeisenbergN9_F1PalindromicPairingViaBlockSpectrum()
    {
        const int N = 9;
        const double J = 1.0;
        const double Gamma = 0.5;
        const double Tolerance = 1e-5;
        double sigma = N * Gamma;  // = 4.5
        const string Topology = "chain (8 bonds)";
        const string JsonFileName = "chain_N9.json";

        // LP64 bridge in place: ComputeSpectrumPerBlock auto-routes blocks larger than
        // Lp64ComplexCeiling (11 585²) through MklDirect + NativeMemory + ILP64-aware LAPACK,
        // bypassing the MathNet marshaller's 2 GB single-array cap. The N=9 chain max block
        // is 15 876² ≈ 4 GB, comfortably above the ceiling so the MklDirect branch carries
        // the load.
        var bonds = ChainBonds(N);

        var totalSw = Stopwatch.StartNew();

        var H = BuildHeisenbergGraphHamiltonian(N, bonds, J: J);
        var gammaPerSite = Enumerable.Repeat(Gamma, N).ToArray();

        var computeSw = Stopwatch.StartNew();
        var spectrum = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N);
        computeSw.Stop();

        Assert.Equal(1 << (2 * N), spectrum.Length);

        var shifted = new Complex[spectrum.Length];
        for (int i = 0; i < spectrum.Length; i++)
            shifted[i] = -2.0 * sigma - spectrum[i];
        MultisetAssert.NearestNeighbourEqual(spectrum, shifted, tolerance: Tolerance,
            context: $"{Topology} N=9 Heisenberg F1 palindromic pairing");

        totalSw.Stop();

        var bondPairs = bonds.Select(b => (b.Site1, b.Site2)).ToArray();
        var metrics = F1SpectrumStatistics.Compute(
            N: N,
            topologyName: Topology,
            hamiltonianClass: "Heisenberg XXX (XX+YY+ZZ) at J=1.0, uniform Z-dephasing γ=0.5",
            jValue: J,
            gammaValue: Gamma,
            bonds: bondPairs,
            spectrum: spectrum,
            totalWallSeconds: totalSw.Elapsed.TotalSeconds,
            computeSpectrumWallSeconds: computeSw.Elapsed.TotalSeconds);

        LogMetrics(metrics, Tolerance);

        string outDir = ResolveMetricsDirectory();
        F1SpectrumStatistics.WriteJson(metrics, outDir, JsonFileName);
        _out.WriteLine($"metrics_json_path={Path.Combine(outDir, JsonFileName)}");
    }

    /// <summary>Mirrors <c>F1GeneralTopologyN8BlockSpectrumTests.LogMetrics</c> — one
    /// grep-able line per metric so a downstream script can parse the test output without
    /// needing the JSON file.</summary>
    private void LogMetrics(F1SpectrumStatistics.TopologyMetrics m, double tolerance)
    {
        _out.WriteLine($"system: N={m.N} topology=\"{m.TopologyName}\" H={m.HamiltonianClass} J={m.JValue} γ={m.GammaValue} σ_shift={m.SigmaShift}");
        _out.WriteLine($"metric=total_wall_seconds value={m.TotalWallSeconds:F3}");
        _out.WriteLine($"metric=compute_spectrum_wall_seconds value={m.ComputeSpectrumWallSeconds:F3}");
        _out.WriteLine($"metric=effective_speedup_over_dense value={m.EffectiveSpeedupOverDense:F1}");
        _out.WriteLine($"metric=tolerance value={tolerance:E0}");
        _out.WriteLine($"metric=max_pairing_distance value={m.MaxPairingDistance:E3}");
        _out.WriteLine($"metric=mean_pairing_distance value={m.MeanPairingDistance:E3}");
        _out.WriteLine($"metric=median_pairing_distance value={m.MedianPairingDistance:E3}");
        _out.WriteLine($"metric=p99_pairing_distance value={m.P99PairingDistance:E3}");
        _out.WriteLine($"metric=min_pairing_distance value={m.MinPairingDistance:E3}");
        _out.WriteLine($"metric=outlier_pair_count value={m.OutlierPairCount}");
        _out.WriteLine($"metric=spectrum_size value={m.SpectrumSize}");
        _out.WriteLine($"metric=min_real value={m.MinReal:E6}");
        _out.WriteLine($"metric=max_real value={m.MaxReal:E6}");
        _out.WriteLine($"metric=min_imag value={m.MinImag:E6}");
        _out.WriteLine($"metric=max_imag value={m.MaxImag:E6}");
        _out.WriteLine($"metric=dissipation_gap value={m.DissipationGap:E6}");
        _out.WriteLine($"metric=kernel_dimension value={m.KernelDimension}");
        _out.WriteLine($"metric=pure_imaginary_count value={m.PureImaginaryCount}");
        _out.WriteLine($"metric=real_eigenvalue_count value={m.RealEigenvalueCount}");
        _out.WriteLine($"metric=distinct_binned_eigenvalue_count value={m.DistinctBinnedEigenvalueCount}");
        _out.WriteLine($"metric=sector_count value={m.SectorCount}");
        _out.WriteLine($"metric=primary_sector_count value={m.PrimarySectorCount}");
        _out.WriteLine($"metric=max_block_size value={m.MaxBlockSize} sector=(p_c={m.MaxBlockSectorPCol},p_r={m.MaxBlockSectorPRow})");
        _out.WriteLine($"metric=top3_block_sizes value=[{string.Join(",", m.Top3BlockSizes)}]");
        _out.WriteLine($"metric=total_block_cubic_cost value={m.TotalBlockCubicCost}");
        _out.WriteLine($"N={m.N} Heisenberg ({m.TopologyName}): {m.SpectrumSize} eigenvalues, " +
                       $"{m.TotalWallSeconds:F1}s total, σ_shift = {m.SigmaShift}, tolerance = {tolerance:E0}, " +
                       $"max pairing distance = {m.MaxPairingDistance:E3}, palindromic pairing OK");
    }

    /// <summary>Same resolution pattern as the N=8 test class: walk up from
    /// <c>AppContext.BaseDirectory</c> to locate <c>simulations/results/</c>, then
    /// descend into the shared <c>f1_n8_n9_metrics/</c> subdirectory.</summary>
    private static string ResolveMetricsDirectory()
    {
        var dir = new DirectoryInfo(AppContext.BaseDirectory);
        while (dir != null)
        {
            var candidate = Path.Combine(dir.FullName, "simulations", "results");
            if (Directory.Exists(candidate))
                return Path.Combine(candidate, "f1_n8_n9_metrics");
            dir = dir.Parent;
        }
        throw new DirectoryNotFoundException(
            $"Cannot locate simulations/results/ by walking up from {AppContext.BaseDirectory}.");
    }
}
