using System.Numerics;

namespace RCPsiSquared.Core.Tests.TestHelpers;

/// <summary>Loader for the EQ-014 precomputed N=7 Liouvillian ground truth produced by
/// <c>compute/RCPsiSquared.Compute</c>'s <c>ptf</c> mode.
///
/// <para>Source: <c>RCPsiSquared.Compute.Program.RunPtfExport</c> — runs the full
/// <c>16384 × 16384</c> dense LAPACK eigendecomposition (left + right eigenvectors,
/// ~16 GB peak memory, ~48 h wall time on commodity hardware) for the chain XY +
/// uniform Z-dephasing Liouvillian at N=7, J=1, γ=0.05. The eigenvalues and biorthogonal
/// eigenvector matrices are written to <c>simulations/results/eq014_*.bin</c> as raw
/// little-endian complex128 in column-major order.</para>
///
/// <para>This helper currently exposes the eigenvalue file (262 KB); the two 4.3 GB
/// eigenvector files are not loaded here — they require streaming access and a different
/// API contract. Eigenvalue cross-validation is enough to verify that any per-sector or
/// reduced-spectrum primitive (e.g. <see cref="RCPsiSquared.Core.BlockSpectrum.LiouvillianBlockSpectrum"/>,
/// <see cref="RCPsiSquared.Core.BlockSpectrum.LiouvillianSectorSweep"/>) reproduces the
/// full N=7 spectrum bit-for-bit.</para>
///
/// <para>Format from <c>eq014_metadata.json</c>:
/// <list type="bullet">
///   <item>N=7, d=128, d²=16384, J=1, γ=0.05, topology=chainXY (PTF convention H = (J/2)(XX+YY))</item>
///   <item>eigenvalue_count = 16384</item>
///   <item>dtype = complex128 little-endian</item>
///   <item>layout = column-major (irrelevant for the 1-D eigenvalue array)</item>
/// </list></para>
/// </summary>
public static class Eq014GroundTruth
{
    public const int N = 7;
    public const int LiouvilleDim = 16384;
    public const double J = 1.0;
    public const double Gamma = 0.05;
    public const string Topology = "chainXY";

    public const string EigenvaluesFileName = "eq014_eigvals_n7.bin";
    private const long EigenvaluesFileSize = (long)LiouvilleDim * 16;  // complex128 = 16 bytes

    /// <summary>Path to <c>simulations/results/</c> derived by walking up from the test
    /// assembly's base directory until that subdirectory exists. Cached on first call.</summary>
    public static string ResultsDirectory => _resultsDir ??= FindResultsDirectory();
    private static string? _resultsDir;

    public static string FullPath(string fileName) => Path.Combine(ResultsDirectory, fileName);

    public static bool IsAvailable(string fileName) =>
        File.Exists(FullPath(fileName));

    /// <summary>Read all <see cref="LiouvilleDim"/> Liouvillian eigenvalues from
    /// <c>eq014_eigvals_n7.bin</c>. Throws <see cref="FileNotFoundException"/> if the file
    /// is absent; tests that depend on it should guard with <see cref="IsAvailable"/> and
    /// skip gracefully when the precomputed file has not been generated locally.</summary>
    public static Complex[] LoadEigenvalues()
    {
        string path = FullPath(EigenvaluesFileName);
        if (!File.Exists(path))
            throw new FileNotFoundException(
                $"Eq014 N=7 eigenvalue ground-truth file missing. " +
                $"Generate via `cd compute/RCPsiSquared.Compute && dotnet run -c Release -- ptf`.",
                path);
        var fi = new FileInfo(path);
        if (fi.Length != EigenvaluesFileSize)
            throw new InvalidDataException(
                $"Eq014 N=7 eigenvalue file size {fi.Length} != expected {EigenvaluesFileSize} " +
                $"({LiouvilleDim} complex128 entries). File may be corrupt or partially written.");

        var eigs = new Complex[LiouvilleDim];
        using var stream = File.OpenRead(path);
        using var reader = new BinaryReader(stream);
        for (int i = 0; i < LiouvilleDim; i++)
        {
            double re = reader.ReadDouble();
            double im = reader.ReadDouble();
            eigs[i] = new Complex(re, im);
        }
        return eigs;
    }

    private static string FindResultsDirectory()
    {
        var dir = new DirectoryInfo(AppContext.BaseDirectory);
        while (dir != null)
        {
            var candidate = Path.Combine(dir.FullName, "simulations", "results");
            if (Directory.Exists(candidate)) return candidate;
            dir = dir.Parent;
        }
        throw new DirectoryNotFoundException(
            $"Cannot locate simulations/results/ by walking up from {AppContext.BaseDirectory}.");
    }
}
