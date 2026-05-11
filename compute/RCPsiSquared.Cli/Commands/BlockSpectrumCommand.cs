using System.Diagnostics;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Cli.Commands;

/// <summary>CLI demonstrator for the BlockSpectrum infrastructure.
///
/// <para>Builds a chain XY+Z-dephasing Liouvillian L at user-specified N, prints sector
/// summary, computes the spectrum either via the full L or via per-block eigendecomposition
/// over the joint-popcount sectors (optionally further refined by the F71 spatial-mirror
/// Z₂), and reports timing + a representative slice of the spectrum.</para>
///
/// <para>Usage: <c>rcpsi block-spectrum --N 6 [--gamma 0.5] [--J 1.0] [--refine f71|none] [--verify]</c>.</para>
///
/// <para>The <c>--verify</c> flag additionally diagonalises the full L directly and asserts
/// that the per-block spectrum matches as a multiset to within 1e-9.</para>
///
/// <para>The F1 palindrome check (set-equality of the spectrum under <c>λ → −2Σγ − λ</c>)
/// is always emitted as the final structural witness.</para></summary>
public static class BlockSpectrumCommand
{
    public static int Run(string[] args)
    {
        var p = new ArgParser(args);
        p.RequireNoPositional();
        int N = p.RequireInt("N");
        double gamma = p.OptionalDouble("gamma") ?? 0.5;
        double J = p.OptionalDouble("J") ?? 1.0;
        string refineKind = p.OptionalString("refine") ?? "none";
        bool verify = p.HasFlag("verify");

        if (N < 1 || N > 12)
            throw new ArgumentOutOfRangeException(nameof(N), N, "Supported N range: 1..12.");
        if (refineKind != "none" && refineKind != "f71")
            throw new ArgumentException($"unknown --refine value: {refineKind}; expected 'none' or 'f71'.");

        Console.WriteLine($"# block-spectrum: N={N}, J={J:G6}, gamma={gamma:G6}, refine={refineKind}");

        // Sector summary (refinement-independent)
        int sectorCount = JointPopcountSectors.SectorCount(N);
        long maxSectorSize = JointPopcountSectors.MaxSectorSize(N);
        long maxBlockMemBytes = maxSectorSize * maxSectorSize * 16L; // complex doubles
        long fullDim = 1L << (2 * N);
        long fullMemBytes = fullDim * fullDim * 16L;

        // Cubic-cost speedup: (4^N)^3 / Σ block_size^3
        double fullCubic = Math.Pow(fullDim, 3);
        double blockCubic = 0;
        for (int pc = 0; pc <= N; pc++)
            for (int pr = 0; pr <= N; pr++)
            {
                long s = JointPopcountSectors.SectorSize(N, pc, pr);
                blockCubic += Math.Pow(s, 3);
            }
        double speedup = fullCubic / blockCubic;

        Console.WriteLine($"# joint-popcount sectors: {sectorCount} = (N+1)^2");
        Console.WriteLine($"# max-block size: {maxSectorSize} (full dim 4^N = {fullDim})");
        Console.WriteLine($"# max-block memory: {FormatBytes(maxBlockMemBytes)} (full L: {FormatBytes(fullMemBytes)})");
        Console.WriteLine($"# cubic-cost speedup vs full eig: ~{speedup:F1}x");

        // Build chain XY + Z-dephasing L
        var sw = Stopwatch.StartNew();
        var L = BuildXYZDephasingL(N, J, gamma);
        sw.Stop();
        Console.WriteLine($"# built L in {sw.ElapsedMilliseconds} ms");

        // Compute spectrum via the chosen path
        Complex[] spectrum;
        sw.Restart();
        if (refineKind == "f71")
        {
            var baseDecomp = JointPopcountSectorBuilder.Build(N);
            var refined = F71MirrorBlockRefinement.RefineWithF71(baseDecomp);
            int maxSubBlock = refined.SectorRanges.Max(s => s.Size);
            int nonEmpty = refined.SectorRanges.Count(s => s.Size > 0);
            Console.WriteLine($"# F71 refinement: {nonEmpty} non-empty sub-blocks (of {refined.SectorRanges.Count} total entries)");
            Console.WriteLine($"# F71-refined max sub-block size: {maxSubBlock}");
            spectrum = F71MirrorBlockRefinement.ComputeSpectrum(L, N);
        }
        else
        {
            spectrum = LiouvillianBlockSpectrum.ComputeSpectrum(L, N);
        }
        sw.Stop();
        Console.WriteLine($"# computed spectrum in {sw.ElapsedMilliseconds} ms ({spectrum.Length} eigenvalues)");

        // Optional verification: full-L direct eig
        if (verify)
        {
            sw.Restart();
            var eigsFull = L.Evd().EigenValues.ToArray();
            sw.Stop();
            Console.WriteLine($"# full-L direct eig in {sw.ElapsedMilliseconds} ms ({eigsFull.Length} eigenvalues)");
            bool ok = MultisetMatches(eigsFull, spectrum, tol: 1e-9, out double maxDev);
            Console.WriteLine($"# verify: multiset match = {(ok ? "PASS" : "FAIL")} (max |Δλ| = {maxDev:E3}, tol = 1e-9)");
            if (!ok) return 3;
        }

        // F1 palindrome check: spectrum invariant under λ → −2·Σγ − λ
        // Σγ = N · gamma (uniform per-site Z-dephasing).
        double sumGamma = N * gamma;
        bool palindrome = MultisetPalindromeOk(spectrum, sumGamma, tol: 1e-7, out double palMaxDev);
        Console.WriteLine($"# F1 palindrome check (λ → −2Σγ − λ; Σγ = {sumGamma:G6}): {(palindrome ? "PASS" : "FAIL")} (max |Δλ| = {palMaxDev:E3})");

        // First ~20 eigenvalues sorted by Re
        Console.WriteLine("# first 20 eigenvalues (sorted by Re ascending):");
        var sorted = spectrum.OrderBy(z => z.Real).ThenBy(z => z.Imaginary).Take(20).ToArray();
        Console.WriteLine("#   idx        Re                Im");
        for (int i = 0; i < sorted.Length; i++)
            Console.WriteLine($"  {i,4}  {sorted[i].Real,17:E9}  {sorted[i].Imaginary,17:E9}");

        return 0;
    }

    private static ComplexMatrix BuildXYZDephasingL(int N, double J, double gamma)
    {
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();
        return PauliDephasingDissipator.BuildZ(H, gammaPerSite);
    }

    private static string FormatBytes(long bytes)
    {
        if (bytes < 1024L) return $"{bytes} B";
        if (bytes < 1024L * 1024) return $"{bytes / 1024.0:F2} KB";
        if (bytes < 1024L * 1024 * 1024) return $"{bytes / (1024.0 * 1024):F2} MB";
        return $"{bytes / (1024.0 * 1024 * 1024):F2} GB";
    }

    /// <summary>Greedy nearest-neighbour multiset matching, identical pattern to the
    /// BlockSpectrum tests' AssertMultisetEqual.</summary>
    private static bool MultisetMatches(IReadOnlyList<Complex> expected, IReadOnlyList<Complex> actual, double tol, out double maxDev)
    {
        maxDev = 0;
        if (expected.Count != actual.Count) return false;
        int n = expected.Count;
        var taken = new bool[n];
        for (int i = 0; i < n; i++)
        {
            double bestDist = double.PositiveInfinity;
            int bestJ = -1;
            for (int j = 0; j < n; j++)
            {
                if (taken[j]) continue;
                double dist = (expected[i] - actual[j]).Magnitude;
                if (dist < bestDist) { bestDist = dist; bestJ = j; }
            }
            if (bestJ < 0) return false;
            if (bestDist > maxDev) maxDev = bestDist;
            if (bestDist >= tol) return false;
            taken[bestJ] = true;
        }
        return true;
    }

    /// <summary>F1 palindromic-spectrum check: for each λ in the spectrum, verify that
    /// its mirror image −2·Σγ − λ also appears in the spectrum (within tol). Greedy
    /// nearest-neighbour matching, mirror image computed on Re only (Im invariant under
    /// reflection of Re).</summary>
    private static bool MultisetPalindromeOk(IReadOnlyList<Complex> spectrum, double sumGamma, double tol, out double maxDev)
    {
        maxDev = 0;
        int n = spectrum.Count;
        var mirrored = new Complex[n];
        for (int i = 0; i < n; i++)
            mirrored[i] = new Complex(-2 * sumGamma - spectrum[i].Real, spectrum[i].Imaginary);
        return MultisetMatches(spectrum, mirrored, tol, out maxDev);
    }
}
