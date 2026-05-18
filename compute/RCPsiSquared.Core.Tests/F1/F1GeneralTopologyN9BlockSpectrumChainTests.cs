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
/// <para><b>Infrastructure ceiling discovered 2026-05-18.</b> At N=9 the max joint-popcount
/// block is C(9, 4) · C(9, 5) = 15 876, i.e. 15 876² × 16 B ≈ 4 GB. MathNet's LP64
/// <c>MklLinearAlgebraProvider.EigenDecomp</c> P/Invoke calls
/// <c>System.StubHelpers.MngdNativeArrayMarshaler.ConvertSpaceToNative</c>, which enforces
/// a hard 2 GB single-native-array limit independent of the <c>AllowVeryLargeObjects</c>
/// CLR flag. The block-spectrum path therefore caps at N=8 on the LP64 MKL route. A
/// 2026-05-18 attempt to run this test threw a paired
/// <c>System.ArgumentException : Array size exceeds addressing limitations</c> from the
/// marshaller after ~1m 44s of setup and the first eig attempt on the largest sector.</para>
///
/// <para><b>This test is marked <see cref="SkippableFactAttribute"/> rather than removed,</b>
/// so it stays discoverable for the future ILP64 bridge. When N=9 becomes reachable
/// (e.g. via <c>RCPsiSquared.Compute.MklDirect</c>'s NativeMemory + ILP64 LAPACK route
/// promoted into Core or invoked from a parallel test fixture), the
/// <see cref="Skip.If"/> call below can be lifted. The test class is kept under
/// <c>[Trait("Category", "SLOW_N9")]</c> so it remains opt-in even after the ceiling is
/// bridged.</para>
///
/// <para>The infrastructure ceiling is the load-bearing finding from the 2026-05-18 SLOW_N8
/// + SLOW_N9 sweep: the F1 palindromic-pairing identity itself is fine at every finite N;
/// what changed at N=9 is that the LP64 MKL P/Invoke marshaller cannot host the largest
/// per-block Evd call. See [PROOF_F1_GENERAL_TOPOLOGY.md § Scale frontier] for the
/// proof-side write-up.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md</c> (scale-frontier section
/// updated with the infrastructure-ceiling finding),
/// <see cref="F1GeneralTopologyVerifiedClaim"/> (typed claim's <c>VerifiedNValues</c>
/// stays at {5, 6, 7, 8}; the <c>ScaleFrontierBlockedAt</c> field documents the N=9
/// blocker), <see cref="F1SpectrumStatistics"/> (shared metrics utility, ready to
/// receive the N=9 capture once the bridge lands),
/// <c>compute/RCPsiSquared.Compute/MklDirect.cs</c> (the NativeMemory + ILP64 LAPACK
/// route that bridges past LP64 marshalling, used by
/// <c>RCPsiSquared.Compute.Liouvillian.BuildDirectNative</c> at N=8 dense and would
/// extend to N=9 per-block).</para></summary>
public class F1GeneralTopologyN9BlockSpectrumChainTests
{
    private readonly ITestOutputHelper _out;

    public F1GeneralTopologyN9BlockSpectrumChainTests(ITestOutputHelper output) => _out = output;

    /// <summary>Hard ceiling on per-block Evd dimension for MathNet's LP64 MKL provider.
    /// The native marshaller <c>MngdNativeArrayMarshaler.ConvertSpaceToNative</c> caps
    /// individual <c>Complex[]</c> P/Invoke arrays at 2 GB; with 16 B per Complex that is
    /// 134 217 728 elements, i.e. a 11 585² square matrix. Block size 15 876² (N=9 max)
    /// exceeds this; block size 4 900² (N=8 max) is comfortably below.</summary>
    public const int Lp64EvdSquareMatrixCeiling = 11_585;

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

    [SkippableFact]
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

        // Skip pre-flight: the LP64 MKL P/Invoke marshaller caps at 11 585² per Evd call.
        // At N=9 the max joint-popcount block is 15 876², exceeding the ceiling. Skipping
        // surfaces the infrastructure gap without aborting the test run; the test stays
        // discoverable so an ILP64 bridge can re-enable it.
        int maxBlockSize = (int)JointPopcountSectors.MaxSectorSize(N);
        long maxBlockBytes = (long)maxBlockSize * maxBlockSize * 16;
        double maxBlockGb = maxBlockBytes / (double)(1L << 30);
        Skip.If(maxBlockSize > Lp64EvdSquareMatrixCeiling,
            $"N={N} max joint-popcount block is {maxBlockSize}² = {(long)maxBlockSize * maxBlockSize:N0} " +
            $"complex elements ({maxBlockGb:F2} GB), exceeding the LP64 MKL P/Invoke 2 GB single-array " +
            $"marshalling ceiling ({Lp64EvdSquareMatrixCeiling}² max square matrix). 2026-05-18 finding: " +
            $"the block-spectrum path caps at N=8 on the LP64 route; reaching N=9 requires routing the " +
            $"dominant block through RCPsiSquared.Compute.MklDirect's NativeMemory + ILP64 LAPACK path " +
            $"(see compute/RCPsiSquared.Compute/MklDirect.cs). Until that bridge lands this test is " +
            $"skipped, not failed: the F1 palindromic-pairing identity itself is exact at every finite N; " +
            $"only the numerical machinery hits the ceiling here.");

        // The remainder of the test runs only once the ceiling is bridged. Implementation is
        // intentionally complete so a future enabling of the path triggers immediate full
        // metric capture without further code change.
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
