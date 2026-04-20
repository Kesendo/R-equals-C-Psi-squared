using System.Diagnostics;
using System.Numerics;
using System.Runtime.InteropServices;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Compute;

// Usage: dotnet run -c Release              -> full suite (N=2-7 + topology + stress + N=8)
//        dotnet run -c Release -- n8        -> N=8 only (skip everything else)
//        dotnet run -c Release -- validate  -> N=5 eigenvalue-only validation only
//        dotnet run -c Release -- rmt       -> RMT eigenvalue export (N=2-7, CSV)
//        dotnet run -c Release -- eigvec    -> Eigenvector export + Pauli projection (N=2-6)
//        dotnet run -c Release -- ptf       -> PTF full dense eigendecomp at N=7 (values + left + right)
bool cavityMode = args.Any(a => a.Equals("cavity", StringComparison.OrdinalIgnoreCase));
bool n8Only = args.Any(a => a.Equals("n8", StringComparison.OrdinalIgnoreCase));
bool validateOnly = args.Any(a => a.Equals("validate", StringComparison.OrdinalIgnoreCase));
bool rmtMode = args.Any(a => a.Equals("rmt", StringComparison.OrdinalIgnoreCase));
bool eigvecMode = args.Any(a => a.Equals("eigvec", StringComparison.OrdinalIgnoreCase));
bool lensMode = args.Any(a => a.Equals("lens", StringComparison.OrdinalIgnoreCase));
bool ptfMode = args.Any(a => a.Equals("ptf", StringComparison.OrdinalIgnoreCase));

bool mklAvailable = false;
try
{
    Control.UseNativeMKL();
    mklAvailable = true;
    Console.WriteLine("MKL Provider: ACTIVE (Intel native LAPACK)");
}
catch (Exception ex)
{
    Console.WriteLine($"MKL Provider: not available ({ex.GetType().Name}), using managed fallback");
}

var resultsDir = Path.Combine(
    @"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared",
    "simulations", "results");

if (cavityMode)
{
    bool testsOnly = args.Any(a => a.Equals("tests", StringComparison.OrdinalIgnoreCase));
    if (testsOnly)
        RunCavityTests(resultsDir);
    else
        RunCavityModes(resultsDir, mklAvailable);
    return;
}

if (rmtMode)
{
    // Optional topology argument: rmt chain (default), rmt star, rmt ring, rmt complete
    string topoArg = args.Where(a => !a.Equals("rmt", StringComparison.OrdinalIgnoreCase))
                         .FirstOrDefault() ?? "chain";
    RunRmtExport(resultsDir, mklAvailable, topoArg.ToLowerInvariant());
    return;
}

if (eigvecMode)
{
    RunEigvecExport(resultsDir);
    return;
}

if (lensMode)
{
    RunLensSurvey(resultsDir);
    return;
}

if (ptfMode)
{
    RunPtfExport(resultsDir);
    return;
}

var finalPath = Path.Combine(resultsDir, "csharp_compute.txt");
var tempPath = Path.Combine(resultsDir, $"csharp_compute_{DateTime.Now:yyyyMMdd_HHmmss}.txt");

using var writer = new StreamWriter(tempPath);

void Log(string msg)
{
    Console.WriteLine(msg);
    writer.WriteLine(msg);
    writer.Flush();
}

const double gamma = 0.05;
var sw = new Stopwatch();
var rateCounts = new Dictionary<int, int>();

Log("=" + new string('=', 79));
Log("R=CPsi^2 COMPUTE ENGINE (C#/.NET + MKL)");
Log($"Started: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
Log($"Machine: {Environment.ProcessorCount} cores, {GC.GetGCMemoryInfo().TotalAvailableMemoryBytes / 1e9:F1} GB RAM");
Log("=" + new string('=', 79));

if (!n8Only && !validateOnly)
{
// ============================================================
// BENCHMARK N=2 to N=7
// ============================================================
Log("\n### BENCHMARK: Liouvillian eigendecomposition timing");
Log($"{"N",4} | {"Matrix",8} | {"Build(ms)",10} | {"Eigen(ms)",10} | {"Total(ms)",10} | {"Rates",7} | {"Mirror",7}");
Log(new string('-', 70));

for (int nQ = 2; nQ <= 7; nQ++)
{
    int d = 1 << nQ;
    int d2 = d * d;
    var bonds = Topology.Star(nQ, Enumerable.Repeat(1.0, nQ - 1).ToArray());
    var gammas = Enumerable.Repeat(gamma, nQ).ToArray();

    List<double> rates;

    if (nQ <= 6)
    {
        // Standard path: MathNet KroneckerProduct + MKL Evd
        sw.Restart();
        var L = Liouvillian.Build(nQ, bonds, gammas);
        var buildMs = sw.ElapsedMilliseconds;

        sw.Restart();
        rates = Liouvillian.GetOscillatoryRates(L);
        var eigenMs = sw.ElapsedMilliseconds;

        var mirror = MirrorAnalysis.CheckSymmetry(rates, nQ * gamma);
        rateCounts[nQ] = rates.Count;
        Log($"{nQ,4} | {d2,8} | {buildMs,10} | {eigenMs,10} | {buildMs + eigenMs,10} | {rates.Count,7} | {mirror.Score,7:P0}");
    }
    else
    {
        // N>=7: Direct element-wise build (no Kronecker) + direct MKL eigen
        Log($"  (N={nQ}: direct build + direct MKL eigen - no Kronecker, no 2GB limit)");

        sw.Restart();
        var rawData = Liouvillian.BuildDirectRaw(nQ, bonds, gammas, msg => Log($"    {msg}"));
        var buildMs = sw.ElapsedMilliseconds;
        Log($"  Build: {buildMs}ms ({buildMs / 1000.0:F1}s)");

        if (mklAvailable)
        {
            Log($"  Eigen via direct MKL (all {Environment.ProcessorCount} cores)...");
            sw.Restart();
            rates = Liouvillian.GetOscillatoryRatesMklRaw(rawData, d2);
            var eigenMs = sw.ElapsedMilliseconds;

            var mirror = MirrorAnalysis.CheckSymmetry(rates, nQ * gamma);
            rateCounts[nQ] = rates.Count;
            Log($"{nQ,4} | {d2,8} | {buildMs,10} | {eigenMs,10} | {buildMs + eigenMs,10} | {rates.Count,7} | {mirror.Score,7:P0}");
        }
        else
        {
            Log($"  MKL not available - skipping N={nQ} (managed too slow)");
        }
    }
}

} // end if (!n8Only && !validateOnly) -- skip N=2-7 benchmark

// ============================================================
// VALIDATION: eigenvalue-only vs z_eigen at N=5
// ============================================================
if (!n8Only && mklAvailable)
{
    Log("\n### VALIDATION: Eigenvalue-only LAPACK path (N=5)");
    int nQv = 5;
    int dv = 1 << nQv;
    int d2v = dv * dv;
    var bondsV = Topology.Star(nQv, Enumerable.Repeat(1.0, nQv - 1).ToArray());
    var gammasV = Enumerable.Repeat(gamma, nQv).ToArray();

    var rawV = Liouvillian.BuildDirectRaw(nQv, bondsV, gammasV);

    // Copy for eigenvalue-only (zgeev destroys input)
    var rawCopy = new Complex[rawV.Length];
    Array.Copy(rawV, rawCopy, rawV.Length);

    // Old path: z_eigen (with eigenvectors)
    sw.Restart();
    var ratesOld = Liouvillian.GetOscillatoryRatesMklRaw(rawV, d2v);
    var msOld = sw.ElapsedMilliseconds;

    // New path: zgeev eigenvalue-only (will use OpenBLAS if MKL zgeev_ unavailable)
    try
    {
        sw.Restart();
        var ratesNew = Liouvillian.GetOscillatoryRatesEigenvaluesOnly(rawCopy, d2v, msg => Log($"  {msg}"));
        var msNew = sw.ElapsedMilliseconds;

        // Compare
        bool match = ratesOld.Count == ratesNew.Count;
        if (match)
        {
            for (int i = 0; i < ratesOld.Count; i++)
            {
                if (Math.Abs(ratesOld[i] - ratesNew[i]) > 0.001)
                { match = false; break; }
            }
        }
        Log($"  z_eigen:   {ratesOld.Count} rates in {msOld}ms");
        Log($"  zgeev(N):  {ratesNew.Count} rates in {msNew}ms");
        Log($"  Match: {(match ? "YES - eigenvalue-only path validated!" : "NO - MISMATCH!")}");

        if (!match)
        {
            Log($"  WARNING: Results differ. Old={ratesOld.Count} New={ratesNew.Count}");
            Log($"  First 5 old: {string.Join(", ", ratesOld.Take(5).Select(r => r.ToString("F6")))}");
            Log($"  First 5 new: {string.Join(", ", ratesNew.Take(5).Select(r => r.ToString("F6")))}");
        }
    }
    catch (Exception ex)
    {
        Log($"  Eigenvalue-only path FAILED: {ex.GetType().Name}: {ex.Message}");
        Log($"  Tried: libMathNetNumericsMKL, mkl_rt.2, libopenblas");
    }
}

// ============================================================
// VALIDATION: ILP64 smoke test at N=5 via native memory path
// ============================================================
if (!n8Only && mklAvailable)
{
    Log("\n### VALIDATION: ILP64 (64-bit int) path via native memory (N=5)");
    int nQi = 5;
    int di = 1 << nQi;
    int d2i = di * di;
    var bondsI = Topology.Star(nQi, Enumerable.Repeat(1.0, nQi - 1).ToArray());
    var gammasI = Enumerable.Repeat(gamma, nQi).ToArray();

    try
    {
        // Build into native memory (same path N=8 will use)
        var ptrI = Liouvillian.BuildDirectNative(nQi, bondsI, gammasI, msg => Log($"  {msg}"));

        // Force ILP64 by calling EigenvaluesOnlyNativeIlp64 through the public API
        // We temporarily lower the threshold - or just call native path directly
        sw.Restart();
        var ratesIlp = Liouvillian.GetOscillatoryRatesNativeIlp64(ptrI, d2i, msg => Log($"  {msg}"));
        var msIlp = sw.ElapsedMilliseconds;

        unsafe { NativeMemory.Free((void*)ptrI); }

        // Compare with known-good LP64 result
        var rawRef = Liouvillian.BuildDirectRaw(nQi, bondsI, gammasI);
        var ratesRef = Liouvillian.GetOscillatoryRatesMklRaw(rawRef, d2i);

        bool match = ratesRef.Count == ratesIlp.Count;
        if (match)
        {
            for (int i = 0; i < ratesRef.Count; i++)
            {
                if (Math.Abs(ratesRef[i] - ratesIlp[i]) > 0.001)
                { match = false; break; }
            }
        }
        Log($"  z_eigen (LP64):       {ratesRef.Count} rates");
        Log($"  zgeev ILP64 (native): {ratesIlp.Count} rates in {msIlp}ms");
        Log($"  Match: {(match ? "YES - ILP64 path validated! N=8 is safe." : "NO - MISMATCH! Do NOT attempt N=8.")}");
    }
    catch (Exception ex)
    {
        Log($"  ILP64 path FAILED: {ex.GetType().Name}: {ex.Message}");
        Log($"  N=8 cannot proceed without a working ILP64 LAPACK.");
    }
}

if (validateOnly)
{
    Log($"\nValidation complete: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
    writer.Close();
    return;
}

// ============================================================
// N=8 ATTEMPT: Native memory + eigenvalue-only LAPACK
// ============================================================
if (mklAvailable)
{
    Log("\n### N=8 ATTEMPT: 65536x65536 matrix (~64 GB)");
    Log($"Available RAM: {GC.GetGCMemoryInfo().TotalAvailableMemoryBytes / 1e9:F1} GB");

    int nQ8 = 8;
    int d8 = 1 << nQ8;
    int d2_8 = d8 * d8;
    long totalElements = (long)d2_8 * d2_8;
    Log($"Matrix: {d2_8}x{d2_8} = {totalElements:N0} elements = {totalElements * 16.0 / 1e9:F1} GB");

    var bonds8 = Topology.Star(nQ8, Enumerable.Repeat(1.0, nQ8 - 1).ToArray());
    var gammas8 = Enumerable.Repeat(gamma, nQ8).ToArray();

    // Ensure OpenBLAS uses all cores for the long eigen
    MklDirect.ConfigureThreads(Environment.ProcessorCount, msg => Log($"  {msg}"));

    try
    {
        // Phase 1: Build into native memory (parallel fill + page warm-up)
        sw.Restart();
        var matrixPtr = Liouvillian.BuildDirectNative(nQ8, bonds8, gammas8, msg => Log($"  {msg}"));
        var buildMs = sw.ElapsedMilliseconds;
        Log($"  Build: {buildMs}ms ({buildMs / 1000.0:F1}s)");

        // Phase 2: Eigenvalues only (no eigenvectors, ILP64 for n>46340)
        Log($"  Eigen via direct LAPACK zgeev ILP64 (eigenvalues only, {Environment.ProcessorCount} cores)...");
        sw.Restart();
        var rates = Liouvillian.GetOscillatoryRatesNative(matrixPtr, d2_8, msg => Log($"  {msg}"));
        var eigenMs = sw.ElapsedMilliseconds;

        // Free native memory immediately
        unsafe { System.Runtime.InteropServices.NativeMemory.Free((void*)matrixPtr); }
        Log($"  Native memory freed.");

        if (rates.Count > 0)
        {
            var mirror = MirrorAnalysis.CheckSymmetry(rates, nQ8 * gamma);
            rateCounts[nQ8] = rates.Count;
            Log($"\n  N=8 RESULT:");
            Log($"  {nQ8,4} | {d2_8,8} | {buildMs,10} | {eigenMs,10} | {buildMs + eigenMs,10} | {rates.Count,7} | {mirror.Score,7:P0}");
            Log($"  Min rate: {rates.Min() / gamma:F3}g  Max rate: {rates.Max() / gamma:F3}g  BW: {(rates.Max() - rates.Min()) / gamma:F3}g");
            Log($"  Palindrome: {mirror.Score:P1} ({mirror.Matched} pairs matched)");
        }
    }
    catch (OutOfMemoryException ex)
    {
        Log($"  OOM: {ex.Message}");
        Log($"  N=8 requires ~64 GB contiguous memory. Close all other applications and retry.");
    }
    catch (Exception ex)
    {
        Log($"  FAILED: {ex.GetType().Name}: {ex.Message}");
        if (ex.Message.Contains("zgeev"))
            Log($"  zgeev_ not found. May need Intel oneAPI MKL (mkl_rt.2.dll) or rebuild libMathNetNumericsMKL with LAPACK exports.");
    }
}
else
{
    Log("\n### N=8: SKIPPED (MKL not available)");
}

if (!n8Only && !validateOnly)
{
// ============================================================
// TOPOLOGY SURVEY at N=4,5,6
// ============================================================
Log("\n### TOPOLOGY SURVEY: Star, Chain, Ring, Complete, Tree");

foreach (int nQ in new[] { 4, 5, 6 })
{
    Log($"\n  N={nQ}:");
    Log($"  {"Topo",10} | {"Rates",7} | {"Min/g",7} | {"Max/g",7} | {"BW/g",7} | {"Mirror",7} | {"ms",7}");
    Log("  " + new string('-', 65));

    var gammas = Enumerable.Repeat(gamma, nQ).ToArray();
    var uniformJ = Enumerable.Repeat(1.0, nQ).ToArray();

    var topos = new (string Name, Bond[] Bonds)[]
    {
        ("star", Topology.Star(nQ, uniformJ[..^1])),
        ("chain", Topology.Chain(nQ, uniformJ[..^1])),
        ("ring", Topology.Ring(nQ, uniformJ)),
        ("complete", Topology.Complete(nQ)),
        ("tree", Topology.BinaryTree(nQ)),
    };

    foreach (var (name, bonds) in topos)
    {
        try
        {
            sw.Restart();
            var L = Liouvillian.Build(nQ, bonds, gammas);
            var rates = Liouvillian.GetOscillatoryRates(L);
            var ms = sw.ElapsedMilliseconds;

            if (rates.Count > 0)
            {
                var mirror = MirrorAnalysis.CheckSymmetry(rates, nQ * gamma);
                Log($"  {name,10} | {rates.Count,7} | {rates.Min() / gamma,7:F3} | {rates.Max() / gamma,7:F3} | " +
                    $"{(rates.Max() - rates.Min()) / gamma,7:F3} | {mirror.Score,7:P0} | {ms,7}");
            }
        }
        catch (Exception ex) { Log($"  {name,10} | FAILED: {ex.Message}"); }
    }
}

// ============================================================
// STRESS TEST: Non-uniform J AND gamma
// ============================================================
Log("\n### STRESS TEST: Non-uniform J AND gamma combined");

foreach (int nQ in new[] { 4, 5, 6 })
{
    var configs = new (double[] J, double[] G)[]
    {
        (Enumerable.Range(0, nQ - 1).Select(k => 0.5 + 0.5 * k).ToArray(),
         Enumerable.Range(0, nQ).Select(k => 0.02 * Math.Pow(2, k)).ToArray()),
        (Enumerable.Range(0, nQ - 1).Select(k => 2.0 - 0.3 * k).ToArray(),
         Enumerable.Range(0, nQ).Select(k => 0.01 + 0.04 * k).ToArray()),
    };

    foreach (var (J, G) in configs)
    {
        var bonds = Topology.Star(nQ, J);
        var L = Liouvillian.Build(nQ, bonds, G);
        var rates = Liouvillian.GetOscillatoryRates(L);
        double sumG = G.Sum();

        if (rates.Count > 0)
        {
            var mid = (rates.Min() + rates.Max()) / 2;
            var symSum = MirrorAnalysis.CheckSymmetry(rates, sumG);
            Log($"  N={nQ} J=[{string.Join(",", J.Select(j => j.ToString("F2")))}] " +
                $"g=[{string.Join(",", G.Select(g => g.ToString("F3")))}]");
            Log($"    sum(g)={sumG:F4}, mid={mid:F5}, sym@sum={symSum.Score:P0}, sum==mid={Math.Abs(sumG - mid) < 0.001}");
        }
    }
}

} // end if (!n8Only) -- skip topology + stress tests

// ============================================================
// SUMMARY
// ============================================================
Log("\n" + new string('=', 80));
Log("FINAL SCALING TABLE");
Log(new string('=', 80));
Log($"{"N",4} | {"Matrix",8} | {"Rates",8} | {"Min/g",7} | {"Max/g",7} | {"BW/g",7} | {"rates/4^N",10}");
Log(new string('-', 60));

foreach (var kvp in rateCounts.OrderBy(x => x.Key))
{
    int n = kvp.Key;
    int count = kvp.Value;
    int d2v = (1 << n) * (1 << n);
    Log($"{n,4} | {d2v,8} | {count,8} | {"2.0",7} | {2 * (n - 1),7:F1} | {2 * (n - 2),7:F1} | {(double)count / Math.Pow(4, n),10:F6}");
}

Log($"\nCompleted: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
Log(new string('=', 80));

writer.Close();
File.Copy(tempPath, finalPath, overwrite: true);
Console.WriteLine($"\n>>> Results saved to: {finalPath}");
Console.WriteLine($"    (timestamped copy: {Path.GetFileName(tempPath)})");

// ============================================================
// RMT EXPORT: All complex eigenvalues as CSV (N=2-7)
// ============================================================
static void RunRmtExport(string resultsDir, bool mklAvailable, string topology = "chain")
{
    Directory.CreateDirectory(resultsDir);
    const double gamma = 0.05;
    const double J = 1.0;

    Console.WriteLine("=== RMT EIGENVALUE EXPORT ===");
    Console.WriteLine($"Date: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
    Console.WriteLine($"Parameters: J={J}, gamma={gamma}, topology={topology}");
    Console.WriteLine();

    var sw = new Stopwatch();

    int nMin = topology == "chain" ? 2 : 3; // non-chain topologies start at N=3
    int nMax = topology == "chain" ? 7 : 6; // non-chain up to N=6

    for (int N = nMin; N <= nMax; N++)
    {
        int d = 1 << N;
        int d2 = d * d;

        Bond[] bonds = topology switch
        {
            "chain" => Topology.Chain(N, Enumerable.Repeat(J, N - 1).ToArray()),
            "star" => Topology.Star(N, Enumerable.Repeat(J, N - 1).ToArray()),
            "ring" => Topology.Ring(N, Enumerable.Repeat(J, N).ToArray()),
            "complete" => Topology.Complete(N, J),
            _ => throw new ArgumentException($"Unknown topology: {topology}. Use chain, star, ring, or complete.")
        };

        var gammas = Enumerable.Repeat(gamma, N).ToArray();

        Console.Write($"N={N} ({d2}x{d2}, {topology})... ");
        sw.Restart();

        Complex[] evals;
        if (N <= 6)
        {
            var L = Liouvillian.Build(N, bonds, gammas);
            evals = Liouvillian.GetAllEigenvalues(L);
        }
        else
        {
            if (!mklAvailable)
            {
                Console.WriteLine("SKIPPED (MKL not available)");
                continue;
            }
            var rawData = Liouvillian.BuildDirectRaw(N, bonds, gammas);
            evals = Liouvillian.GetAllEigenvaluesMklRaw(rawData, d2);
        }

        // Chain uses legacy filenames (no prefix) for backwards compatibility
        string csvName = topology == "chain"
            ? $"rmt_eigenvalues_N{N}.csv"
            : $"rmt_eigenvalues_{topology}_N{N}.csv";
        var csvPath = Path.Combine(resultsDir, csvName);

        using (var csvWriter = new StreamWriter(csvPath))
        {
            csvWriter.WriteLine("Re\tIm");
            foreach (var ev in evals)
            {
                csvWriter.Write(ev.Real.ToString("R", System.Globalization.CultureInfo.InvariantCulture));
                csvWriter.Write('\t');
                csvWriter.WriteLine(ev.Imaginary.ToString("R", System.Globalization.CultureInfo.InvariantCulture));
            }
        }

        Console.WriteLine($"{evals.Length} eigenvalues in {sw.ElapsedMilliseconds}ms -> {Path.GetFileName(csvPath)}");
    }

    Console.WriteLine();
    Console.WriteLine("=== RMT EXPORT COMPLETE ===");
    Console.WriteLine($"Files in: {resultsDir}");
}

// ============================================================
// CAVITY MODES: eigenvalue analysis at zero noise (gamma=0)
// ============================================================
static void RunCavityModes(string resultsDir, bool mklAvailable)
{
    Directory.CreateDirectory(resultsDir);
    var outPath = Path.Combine(resultsDir, "cavity_modes_zero_noise.txt");
    using var outFile = new StreamWriter(outPath);

    void CLog(string s = "")
    {
        Console.WriteLine(s);
        outFile.WriteLine(s);
        outFile.Flush();
    }

    CLog("=== CAVITY MODES AT ZERO NOISE ===");
    CLog($"Date: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
    CLog();

    // Python reference values (from ZERO_IS_THE_MIRROR.md)
    var pyRef = new Dictionary<int, (int stat, int osc, int freqs)>
    {
        [2] = (10, 6, 1),
        [3] = (24, 40, 3),
        [4] = (54, 202, 14),
        [5] = (120, 904, 43),
    };

    double J = 1.0;
    var results = new List<(int N, int d2, Liouvillian.CavityModes star,
                            Liouvillian.CavityModes? chain)>();

    foreach (int N in new[] { 2, 3, 4, 5, 6, 7 })
    {
        int d = 1 << N;
        int d2 = d * d;
        var gammas = new double[N]; // all zeros

        // Use Chain for Python verification (ZERO_IS_THE_MIRROR used chain)
        var chainBondsMain = Topology.Chain(N, Enumerable.Repeat(J, N - 1).ToArray());

        CLog($"N={N}: Chain topology, J={J}, gamma=0.0");
        CLog($"  Liouvillian: {d2}x{d2}");

        var sw = Stopwatch.StartNew();
        Liouvillian.CavityModes starModes;

        if (N <= 6)
        {
            var L = Liouvillian.Build(N, chainBondsMain, gammas);
            starModes = Liouvillian.GetCavityModes(L);
        }
        else
        {
            if (!mklAvailable)
            {
                CLog("  MKL not available, skipping N=7");
                CLog();
                continue;
            }
            var rawData = Liouvillian.BuildDirectRaw(N, chainBondsMain, gammas, s => CLog($"  {s}"));
            starModes = Liouvillian.GetCavityModesMklRaw(rawData, d2);
        }

        CLog($"  Total eigenvalues: {starModes.Total}");
        CLog($"  Stationary (Re=0, Im=0): {starModes.Stationary}");
        CLog($"  Oscillating (Re=0, Im!=0): {starModes.Oscillating}");
        if (starModes.Anomalous > 0)
            CLog($"  WARNING: Anomalous (Re!=0): {starModes.Anomalous}");
        CLog($"  Max |Re|: {starModes.MaxAbsReal:E2}");
        CLog($"  Distinct frequencies: {starModes.Frequencies.Length}");

        if (starModes.Frequencies.Length <= 50)
        {
            var freqStr = string.Join(", ", starModes.Frequencies.Select(f => $"{f / J:F3}"));
            CLog($"  Frequencies / J: [{freqStr}]");

            var intMultiples = starModes.Frequencies.Select(f => f / (2.0 * J)).ToArray();
            var allInteger = intMultiples.All(m => Math.Abs(m - Math.Round(m)) < 0.01);
            if (allInteger && starModes.Frequencies.Length > 0)
                CLog($"  Integer multiples of 2J: [{string.Join(", ", intMultiples.Select(m => $"{Math.Round(m)}"))}]");
        }

        // Verify against Python
        if (pyRef.TryGetValue(N, out var expected))
        {
            bool statPass = starModes.Stationary == expected.stat;
            bool oscPass = starModes.Oscillating == expected.osc;
            bool freqPass = starModes.Frequencies.Length == expected.freqs;
            string status = (statPass && oscPass && freqPass) ? "PASS" : "FAIL";
            CLog($"  VERIFY vs Python: {starModes.Stationary}+{starModes.Oscillating} " +
                 $"(expected {expected.stat}+{expected.osc}), freqs={starModes.Frequencies.Length} " +
                 $"(expected {expected.freqs}) = {status}");

            if (status == "FAIL")
            {
                CLog();
                CLog("*** VERIFICATION FAILED. STOPPING. ***");
                return;
            }
        }
        else
        {
            CLog($"  NEW DATA POINT (no Python reference)");
        }

        CLog($"  Time: {sw.Elapsed}");
        CLog();

        // Star topology comparison for N=3 and N=4
        Liouvillian.CavityModes? altModes = null;
        if (N >= 3 && N <= 6)
        {
            var starBonds = Topology.Star(N, Enumerable.Repeat(J, N - 1).ToArray());
            var swS = Stopwatch.StartNew();
            var Ls = Liouvillian.Build(N, starBonds, gammas);
            altModes = Liouvillian.GetCavityModes(Ls);

            CLog($"N={N}: Star topology, J={J}, gamma=0.0");
            CLog($"  Stationary: {altModes.Stationary}, Oscillating: {altModes.Oscillating}");
            CLog($"  Distinct frequencies: {altModes.Frequencies.Length}");
            if (altModes.Frequencies.Length <= 50)
            {
                var freqStr = string.Join(", ", altModes.Frequencies.Select(f => $"{f / J:F3}"));
                CLog($"  Frequencies / J: [{freqStr}]");
            }
            bool same = altModes.Stationary == starModes.Stationary
                     && altModes.Oscillating == starModes.Oscillating
                     && altModes.Frequencies.Length == starModes.Frequencies.Length;
            CLog($"  Same as Chain? {(same ? "YES" : "NO")}");
            CLog($"  Time: {swS.Elapsed}");
            CLog();
        }

        results.Add((N, d2, starModes, altModes));
    }

    // Summary table
    CLog(new string('=', 80));
    CLog("SUMMARY TABLE (Chain topology)");
    CLog(new string('=', 80));
    CLog();
    CLog("| N | d^2 | Stationary | Oscillating | Frequencies | Max|Re| |");
    CLog("|---|-----|-----------|-------------|-------------|---------|");
    foreach (var (N, d2, m, _) in results)
    {
        CLog($"| {N} | {d2} | {m.Stationary} | {m.Oscillating} | {m.Frequencies.Length} | {m.MaxAbsReal:E1} |");
    }
    CLog();

    // Chain comparison
    var chainResults = results.Where(r => r.chain != null).ToList();
    if (chainResults.Count > 0)
    {
        CLog("CHAIN vs STAR COMPARISON");
        CLog();
        foreach (var (N, _, chain, star) in chainResults)
        {
            CLog($"N={N}: Chain({chain.Stationary}+{chain.Oscillating}, {chain.Frequencies.Length} freq) " +
                 $"vs Star({star!.Stationary}+{star.Oscillating}, {star.Frequencies.Length} freq)");
        }
        CLog();
    }

    CLog("=== DONE ===");
    Console.WriteLine($"\n>>> Results saved to: {outPath}");
}

// ============================================================
// CAVITY TESTS: Ring, Complete, non-uniform J
// ============================================================
static void RunCavityTests(string resultsDir)
{
    Directory.CreateDirectory(resultsDir);
    var outPath = Path.Combine(resultsDir, "cavity_modes_tests.txt");
    using var outFile = new StreamWriter(outPath);

    void CLog(string s = "")
    {
        Console.WriteLine(s);
        outFile.WriteLine(s);
        outFile.Flush();
    }

    CLog("=== CAVITY MODES: TOPOLOGY AND COUPLING TESTS ===");
    CLog($"Date: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
    CLog();

    // Expected stationary counts from Clebsch-Gordan formula
    var expectedStat = new Dictionary<int, int>
    {
        [2] = 10, [3] = 24, [4] = 54, [5] = 120, [6] = 260
    };

    double J = 1.0;

    // ---- TEST A: Ring topology ----
    // Ring has C_N spatial symmetry -> additional degeneracies -> MORE stationary modes
    CLog("### TEST A: Ring topology at gamma=0");
    CLog("(Spatial symmetry may create additional degeneracies: Stat >= formula)");
    CLog();

    foreach (int N in new[] { 3, 4, 5, 6 })
    {
        int d = 1 << N;
        var gammas = new double[N];
        var bonds = Topology.Ring(N, Enumerable.Repeat(J, N).ToArray());
        var sw = Stopwatch.StartNew();
        var L = Liouvillian.Build(N, bonds, gammas);
        var modes = Liouvillian.GetCavityModes(L);

        string verify = expectedStat.TryGetValue(N, out var exp)
            ? (modes.Stationary >= exp ? $"OK (>={exp})" : "FAIL")
            : "NEW";

        CLog($"N={N} Ring: Stat={modes.Stationary}, Osc={modes.Oscillating}, Freq={modes.Frequencies.Length} " +
             $"Max|Re|={modes.MaxAbsReal:E1} [{verify}] ({sw.Elapsed})");

        if (modes.Frequencies.Length <= 20)
        {
            var freqStr = string.Join(", ", modes.Frequencies.Select(f => $"{f / J:F3}"));
            CLog($"  Frequencies / J: [{freqStr}]");
        }
    }
    CLog();

    // ---- TEST B: Complete topology ----
    CLog("### TEST B: Complete topology at gamma=0");
    CLog("(S_N symmetry -> maximal degeneracy -> fewest frequencies)");
    CLog();

    foreach (int N in new[] { 3, 4, 5 })
    {
        int d = 1 << N;
        var gammas = new double[N];
        var bonds = Topology.Complete(N);
        var sw = Stopwatch.StartNew();
        var L = Liouvillian.Build(N, bonds, gammas);
        var modes = Liouvillian.GetCavityModes(L);

        string verify = expectedStat.TryGetValue(N, out var exp)
            ? (modes.Stationary >= exp ? $"OK (>={exp})" : "FAIL")
            : "NEW";

        CLog($"N={N} Complete: Stat={modes.Stationary}, Osc={modes.Oscillating}, Freq={modes.Frequencies.Length} " +
             $"Max|Re|={modes.MaxAbsReal:E1} [{verify}] ({sw.Elapsed})");

        if (modes.Frequencies.Length <= 20)
        {
            var freqStr = string.Join(", ", modes.Frequencies.Select(f => $"{f / J:F3}"));
            CLog($"  Frequencies / J: [{freqStr}]");
        }
    }
    CLog();

    // ---- TEST C: Non-uniform J couplings ----
    CLog("### TEST C: Non-uniform J couplings on Chain N=4");
    CLog();

    var jConfigs = new (string name, double[] jvals)[]
    {
        ("uniform", new[] { 1.0, 1.0, 1.0 }),
        ("linear",  new[] { 0.5, 1.0, 1.5 }),
        ("random",  new[] { 0.3, 2.1, 0.7 }),
    };

    foreach (var (name, jvals) in jConfigs)
    {
        int N = 4;
        var gammas = new double[N];
        var bonds = Topology.Chain(N, jvals);
        var sw = Stopwatch.StartNew();
        var L = Liouvillian.Build(N, bonds, gammas);
        var modes = Liouvillian.GetCavityModes(L);

        string verify = modes.Stationary == 54 ? "PASS" : "FAIL";

        CLog($"N=4 Chain J=[{string.Join(",", jvals.Select(j => j.ToString("F1")))}] ({name}): " +
             $"Stat={modes.Stationary}, Osc={modes.Oscillating}, Freq={modes.Frequencies.Length} [{verify}] ({sw.Elapsed})");

        if (modes.Frequencies.Length <= 20)
        {
            var freqStr = string.Join(", ", modes.Frequencies.Select(f => $"{f:F3}"));
            CLog($"  Frequencies: [{freqStr}]");
        }

        if (verify == "FAIL")
        {
            CLog("*** VERIFICATION FAILED. STOPPING. ***");
            return;
        }
    }
    CLog();

    // ---- SUMMARY ----
    CLog(new string('=', 60));
    CLog("ALL TESTS COMPLETED");
    CLog("Clebsch-Gordan formula Stat(N) = Sum_J m(J,N)*(2J+1)^2");
    CLog("is EXACT for Chain (minimal symmetry) and a LOWER BOUND for");
    CLog("higher-symmetry topologies (Ring, Complete, Star).");
    CLog("Non-uniform J on Chain: formula still exact (SU(2) preserved).");
    CLog(new string('=', 60));

    Console.WriteLine($"\n>>> Results saved to: {outPath}");
}

// ============================================================
// EIGENVECTOR EXPORT: Pauli projection at Re = -2γ (N=2-6)
// ============================================================
static void RunEigvecExport(string resultsDir)
{
    Directory.CreateDirectory(resultsDir);
    const double gamma = 0.05;
    const double J = 1.0;
    const double gridSpacing = 2.0 * gamma; // = 0.1, Lindblad convention
    const double tol = 1e-8;

    Console.WriteLine("=== EIGENVECTOR EXPORT + PAULI PROJECTION ===");
    Console.WriteLine($"Date: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
    Console.WriteLine($"Parameters: J={J}, gamma={gamma}, Chain topology");
    Console.WriteLine($"Target: Re = -{gridSpacing} (first non-zero grid point, weight w=1)");
    Console.WriteLine();

    var sw = new Stopwatch();

    for (int N = 2; N <= 6; N++)
    {
        int d = 1 << N;
        int d2 = d * d;
        var bonds = Topology.Chain(N, Enumerable.Repeat(J, N - 1).ToArray());
        var gammas = Enumerable.Repeat(gamma, N).ToArray();

        Console.WriteLine(new string('=', 70));
        Console.WriteLine($"N={N}: {d2}x{d2} Liouvillian, expecting {2 * N} real eigenvalues at Re=-{gridSpacing}");
        Console.WriteLine(new string('=', 70));

        sw.Restart();
        var L = Liouvillian.Build(N, bonds, gammas);
        var (values, vectors) = Liouvillian.GetAllEigenvaluesAndVectors(L);
        Console.WriteLine($"  Eigendecomposition: {sw.ElapsedMilliseconds}ms, {values.Length} eigenvalues");

        // Find purely real eigenvalues at Re = -gridSpacing
        var targetIndices = new List<int>();
        for (int i = 0; i < values.Length; i++)
        {
            if (Math.Abs(values[i].Real + gridSpacing) < tol && Math.Abs(values[i].Imaginary) < tol)
                targetIndices.Add(i);
        }

        Console.WriteLine($"  Found {targetIndices.Count} eigenvalues at Re=-{gridSpacing}, Im=0 (expected {2 * N})");

        if (targetIndices.Count != 2 * N)
            Console.WriteLine($"  WARNING: count mismatch! Expected {2 * N}, got {targetIndices.Count}");

        // Verification: L*v = λ*v for first and last eigenvector
        if (targetIndices.Count > 0)
        {
            foreach (int checkIdx in new[] { targetIndices[0], targetIndices[^1] })
            {
                var lambda = values[checkIdx];
                var v = new Complex[d2];
                for (int i = 0; i < d2; i++)
                    v[i] = vectors[(long)checkIdx * d2 + i];

                // Compute L*v
                var lv = new Complex[d2];
                for (int i = 0; i < d2; i++)
                    for (int j = 0; j < d2; j++)
                        lv[i] += L[i, j] * v[j];

                // Residual ||L*v - λ*v|| / ||v||
                double normV = 0, normRes = 0;
                for (int i = 0; i < d2; i++)
                {
                    normV += (v[i] * Complex.Conjugate(v[i])).Real;
                    var diff = lv[i] - lambda * v[i];
                    normRes += (diff * Complex.Conjugate(diff)).Real;
                }
                double relRes = Math.Sqrt(normRes / normV);
                string status = relRes < 1e-8 ? "PASS" : "FAIL";
                Console.WriteLine($"  Verification L*v=lv for eigvec #{checkIdx}: residual={relRes:E2} [{status}]");
            }
        }

        // Pauli projection for each target eigenvector
        Console.WriteLine($"\n  Pauli projection of {targetIndices.Count} eigenvectors...");
        sw.Restart();

        var csvPath = Path.Combine(resultsDir, $"eigvec_at_minus_gamma_N{N}.csv");
        using var csvWriter = new StreamWriter(csvPath);
        csvWriter.WriteLine("eigvec_idx\teigenvalue_re\teigenvalue_im\tpauli_string\txy_weight\tcoeff_abs\tcoeff_phase");

        // Track dominant Pauli strings across all eigenvectors
        var dominantStrings = new List<(int eigIdx, string label, int weight, double absCoeff)>();

        for (int ei = 0; ei < targetIndices.Count; ei++)
        {
            int evIdx = targetIndices[ei];
            var lambda = values[evIdx];

            // Extract eigenvector
            var v = new Complex[d2];
            for (int i = 0; i < d2; i++)
                v[i] = vectors[(long)evIdx * d2 + i];

            // Project onto Pauli basis
            var coeffs = PauliOps.ProjectOntoPauliBasis(v, N);

            // Find dominant coefficient and write CSV
            double maxAbs = 0;
            string maxLabel = "";
            int maxWeight = 0;

            for (int s = 0; s < coeffs.Length; s++)
            {
                double absC = coeffs[s].Magnitude;
                if (absC > 1e-10)
                {
                    string label = PauliOps.PauliLabel(s, N);
                    int w = PauliOps.XYWeight(s, N);
                    double phase = coeffs[s].Phase;

                    csvWriter.Write(ei.ToString(System.Globalization.CultureInfo.InvariantCulture));
                    csvWriter.Write('\t');
                    csvWriter.Write(lambda.Real.ToString("R", System.Globalization.CultureInfo.InvariantCulture));
                    csvWriter.Write('\t');
                    csvWriter.Write(lambda.Imaginary.ToString("R", System.Globalization.CultureInfo.InvariantCulture));
                    csvWriter.Write('\t');
                    csvWriter.Write(label);
                    csvWriter.Write('\t');
                    csvWriter.Write(w.ToString());
                    csvWriter.Write('\t');
                    csvWriter.Write(absC.ToString("R", System.Globalization.CultureInfo.InvariantCulture));
                    csvWriter.Write('\t');
                    csvWriter.WriteLine(phase.ToString("R", System.Globalization.CultureInfo.InvariantCulture));
                }

                if (absC > maxAbs)
                {
                    maxAbs = absC;
                    maxLabel = PauliOps.PauliLabel(s, N);
                    maxWeight = PauliOps.XYWeight(s, N);
                }
            }

            dominantStrings.Add((ei, maxLabel, maxWeight, maxAbs));
        }

        Console.WriteLine($"  Pauli projection: {sw.ElapsedMilliseconds}ms -> {Path.GetFileName(csvPath)}");

        // Summary: dominant Pauli strings
        Console.WriteLine($"\n  Dominant Pauli strings (2N={2 * N} eigenvectors at Re=-{gridSpacing}):");
        Console.WriteLine($"  {"#",3} {"Dominant",10} {"w",3} {"|c|",10}");
        Console.WriteLine("  " + new string('-', 30));
        foreach (var (idx, label, w, absC) in dominantStrings)
            Console.WriteLine($"  {idx,3} {label,10} {w,3} {absC,10:F6}");

        // Weight sector distribution
        var weightCounts = dominantStrings.GroupBy(x => x.weight)
            .OrderBy(g => g.Key)
            .Select(g => $"w={g.Key}: {g.Count()}")
            .ToArray();
        Console.WriteLine($"\n  Weight distribution: {string.Join(", ", weightCounts)}");

        // Count distinct labels
        var distinctLabels = dominantStrings.Select(x => x.label).Distinct().Count();
        Console.WriteLine($"  Distinct dominant labels: {distinctLabels}");

        // Detailed: for each eigenvector, show top 3 Pauli coefficients
        Console.WriteLine($"\n  Top-3 Pauli coefficients per eigenvector:");
        for (int ei = 0; ei < targetIndices.Count; ei++)
        {
            int evIdx = targetIndices[ei];
            var v = new Complex[d2];
            for (int i = 0; i < d2; i++)
                v[i] = vectors[(long)evIdx * d2 + i];

            var coeffs = PauliOps.ProjectOntoPauliBasis(v, N);

            // Sort by magnitude
            var ranked = Enumerable.Range(0, coeffs.Length)
                .Select(s => (idx: s, label: PauliOps.PauliLabel(s, N),
                              weight: PauliOps.XYWeight(s, N), abs: coeffs[s].Magnitude))
                .Where(x => x.abs > 1e-10)
                .OrderByDescending(x => x.abs)
                .Take(3)
                .ToArray();

            var parts = ranked.Select(r => $"{r.label}(w={r.weight}, |c|={r.abs:F4})");
            Console.WriteLine($"    #{ei}: {string.Join("  ", parts)}");
        }

        Console.WriteLine();
    }

    Console.WriteLine("=== EIGENVECTOR EXPORT COMPLETE ===");
    Console.WriteLine($"Files in: {resultsDir}");
}

// ============================================================
// LENS SURVEY: slow-mode lens extraction across N, topology, gamma
// ============================================================
static void RunLensSurvey(string resultsDir)
{
    var outDir = Path.Combine(resultsDir, "lens_survey");
    Directory.CreateDirectory(outDir);
    var summaryPath = Path.Combine(outDir, "lens_survey_summary.txt");
    var scalingPath = Path.Combine(outDir, "lens_survey_scaling.txt");
    var jsonPath = Path.Combine(outDir, "lens_survey_results.json");

    Console.WriteLine("=" + new string('=', 79));
    Console.WriteLine("LENS SURVEY: Slow-mode lens extraction");
    Console.WriteLine("=" + new string('=', 79));

    // ---- Validation gate: N=5 IBM sacrifice profile ----
    Console.WriteLine("\n--- VALIDATION: N=5 chain, IBM sacrifice profile ---");
    double[] ibmGamma = { 2.33573, 0.09937, 0.05000, 0.07173, 0.05132 };
    var valBonds = Topology.Chain(5, Enumerable.Repeat(1.0, 4).ToArray());
    var valResult = LensAnalysis.RunFullLensPipeline(
        5, valBonds, ibmGamma, "N5_chain_ibm_sacrifice", "chain",
        msg => { Console.WriteLine(msg); Console.Out.Flush(); });

    // Validation checks
    bool valOk = true;

    // The lens state should be extracted from the SE-accessible slow mode (rate ~0.318)
    if (valResult.LensState != null)
    {
        // Find the mode the lens was extracted from (first with SE > 0.01)
        double lensRate = -valResult.SlowModes[0].Eigenvalue.Real;
        double seFrac = valResult.SEFrobRatios[0];

        Console.WriteLine($"\n  Lens mode rate: {lensRate:F4} (expect 0.3181)");
        if (Math.Abs(lensRate - 0.3181) > 0.001) { Console.WriteLine("  ** FAIL: rate mismatch"); valOk = false; }

        Console.WriteLine($"  SE Frob ratio:  {seFrac:F4} (expect 0.9986)");
        if (Math.Abs(seFrac - 0.9986) > 0.005) { Console.WriteLine("  ** FAIL: SE fraction mismatch"); valOk = false; }

        double[] refAmps = { 0.099, 0.239, 0.428, 0.572, 0.651 };
        double cosine = 0;
        for (int i = 0; i < 5; i++)
            cosine += valResult.LensState.Amplitudes[i] * refAmps[i];
        Console.WriteLine($"  Cosine sim:     {cosine:F6} (expect > 0.999)");
        if (cosine < 0.999) { Console.WriteLine("  ** FAIL: amplitude mismatch"); valOk = false; }

        Console.WriteLine($"  |c_slow|:       {valResult.LensState.SlowModeProjection:F4} (expect 0.972)");
        Console.WriteLine($"  Amplitudes:     [{string.Join(", ", valResult.LensState.Amplitudes.Select(a => a.ToString("F3")))}]");
    }
    else { Console.WriteLine("  ** FAIL: no lens state extracted"); valOk = false; }

    // Check the inaccessible second mode (SecondModeRate is Re(lambda), negative)
    double secondRate = Math.Abs(valResult.SecondModeRate);
    Console.WriteLine($"  2nd mode rate:  {secondRate:F4} (expect 0.1674)");
    double secondSE = valResult.SEFrobRatios.Length > 1 ? valResult.SEFrobRatios[1] : double.NaN;
    Console.WriteLine($"  2nd mode SE:    {secondSE:E2} (expect < 1e-10)");

    if (!valOk)
    {
        Console.WriteLine("\n  VALIDATION FAILED. Stopping.");
        return;
    }
    Console.WriteLine("\n  VALIDATION PASSED. Proceeding to sweep.\n");

    // ---- Build configuration matrix ----
    var configs = new List<(int n, string topo, Func<int, double[], Bond[]> bondGen, string profile, Func<int, double[]> gammaGen)>();

    double gammaBase = 0.05;
    double epsilon = 0.001;

    Func<int, double[]> uniformGamma = n =>
        Enumerable.Repeat(gammaBase, n).ToArray();

    Func<int, double[]> edgeSacrifice = n =>
    {
        var g = new double[n];
        g[0] = n * gammaBase - (n - 1) * epsilon;
        for (int i = 1; i < n; i++) g[i] = epsilon;
        return g;
    };

    Func<int, double[]> centerSacrifice = n =>
    {
        int center = n / 2;
        var g = new double[n];
        g[center] = n * gammaBase - (n - 1) * epsilon;
        for (int i = 0; i < n; i++)
            if (i != center) g[i] = epsilon;
        return g;
    };

    Func<int, double[]> moderateAsymmetry = n =>
    {
        var g = new double[n];
        for (int i = 0; i < n; i++)
            g[i] = gammaBase * (1.0 + 0.5 * i / (n - 1));
        return g;
    };

    var profiles = new (string name, Func<int, double[]> gen)[]
    {
        ("uniform", uniformGamma),
        ("edge_sacrifice", edgeSacrifice),
        ("center_sacrifice", centerSacrifice),
        ("moderate_asymmetry", moderateAsymmetry),
    };

    var topos = new (string name, Func<int, double[], Bond[]> gen)[]
    {
        ("chain", (n, c) => Topology.Chain(n, c)),
        ("star", (n, c) => Topology.Star(n, c)),
        ("ring", (n, c) => Topology.Ring(n, c)),
    };

    foreach (var (topoName, bondGen) in topos)
    {
        int maxN = topoName == "chain" ? 7 : 6; // N=7 chain only (needs gcAllowVeryLargeObjects)
        for (int n = 2; n <= maxN; n++)
        {
            if (topoName == "star" && n < 3) continue;
            foreach (var (profName, gammaGen) in profiles)
            {
                var couplings = Enumerable.Repeat(1.0, n * n).ToArray(); // enough for any topology
                configs.Add((n, topoName, bondGen, profName, gammaGen));
            }
        }
    }

    // Also Complete at N=3,4 as sanity
    Func<int, double[], Bond[]> completeBondGen = (nn, c) => Topology.Complete(nn);
    for (int n = 3; n <= 4; n++)
    {
        foreach (var (profName, gammaGen) in profiles)
        {
            configs.Add((n, "complete", completeBondGen, profName, gammaGen));
        }
    }

    Console.WriteLine($"Total configurations: {configs.Count}");

    // ---- Run sweep ----
    var allResults = new List<LensAnalysis.LensSurveyResult>();
    allResults.Add(valResult); // include validation run

    var sw = System.Diagnostics.Stopwatch.StartNew();
    var errorLogPath = Path.Combine(outDir, "lens_survey_errors.log");

    using var summaryWriter = new StreamWriter(summaryPath);
    summaryWriter.WriteLine("Lens Survey Summary");
    summaryWriter.WriteLine($"Computed: {DateTime.Now:yyyy-MM-dd HH:mm}");
    summaryWriter.WriteLine();
    summaryWriter.WriteLine($"{"N",2}  {"Topology",-10}  {"Profile",-20}  {"SlowRate",10}  {"SEFrac",8}  {"psi_opt_shape",-40}  {"|c_slow|",8}  {"2ndRate",10}  {"2ndSE",10}");
    summaryWriter.WriteLine(new string('-', 130));

    for (int ci = 0; ci < configs.Count; ci++)
    {
        var (n, topoName, bondGen, profName, gammaGen) = configs[ci];
        var gammas = gammaGen(n);
        var nBonds = topoName == "ring" ? n : (topoName == "complete" ? n * (n - 1) / 2 : n - 1);
        var couplings = Enumerable.Repeat(1.0, nBonds).ToArray();

        Bond[] bonds;
        try { bonds = bondGen(n, couplings); }
        catch { Console.WriteLine($"  Skip: {topoName} N={n} (bond generation failed)"); continue; }

        string label = $"N{n}_{topoName}_{profName}";
        Console.Write($"  [{ci + 1}/{configs.Count}] {label,-40}");
        Console.Out.Flush();

        var runSw = System.Diagnostics.Stopwatch.StartNew();
        LensAnalysis.LensSurveyResult result;
        try
        {
            result = LensAnalysis.RunFullLensPipeline(n, bonds, gammas, label, topoName,
                msg => { Console.WriteLine(msg); Console.Out.Flush(); });
        }
        catch (Exception ex)
        {
            Console.WriteLine($" ERROR: {ex.Message}");
            Console.Out.Flush();
            File.AppendAllText(errorLogPath,
                $"[{DateTime.Now:HH:mm:ss}] {label}: {ex.GetType().Name}: {ex.Message}\n{ex.StackTrace}\n\n");
            continue;
        }
        runSw.Stop();
        allResults.Add(result);

        double rate = result.SlowModes.Length > 0 ? -result.SlowModes[0].Eigenvalue.Real : double.NaN;
        double se = result.SEFrobRatios.Length > 0 ? result.SEFrobRatios[0] : double.NaN;
        string shape = result.LensState != null
            ? "[" + string.Join(",", result.LensState.Amplitudes.Select(a => a.ToString("F3"))) + "]"
            : "n/a";
        double cslow = result.LensState?.SlowModeProjection ?? double.NaN;
        double r2 = result.SecondModeRate;
        double se2 = result.SEFrobRatios.Length > 1 ? result.SEFrobRatios[1] : double.NaN;

        Console.WriteLine($" rate={rate:F4} SE={se:F4} |c|={cslow:F3} ({runSw.Elapsed.TotalSeconds:F1}s)");
        Console.Out.Flush();

        summaryWriter.WriteLine($"{n,2}  {topoName,-10}  {profName,-20}  {rate,10:F4}  {se,8:F4}  {shape,-40}  {cslow,8:F4}  {r2,10:F4}  {se2,10:E2}");
        summaryWriter.Flush();
        summaryWriter.Flush();
    }

    sw.Stop();

    // ---- Scaling analysis ----
    using var scalingWriter = new StreamWriter(scalingPath);
    scalingWriter.WriteLine("Lens Survey Scaling Analysis");
    scalingWriter.WriteLine($"Computed: {DateTime.Now:yyyy-MM-dd HH:mm}");

    scalingWriter.WriteLine("\n=== SE FRACTION vs N (edge sacrifice, chain) ===");
    foreach (var r in allResults.Where(r => r.Topology == "chain" && r.Label.Contains("edge_sacrifice")))
        scalingWriter.WriteLine($"N={r.N}: {(r.SEFrobRatios.Length > 0 ? r.SEFrobRatios[0].ToString("F6") : "n/a")}");

    scalingWriter.WriteLine("\n=== MONOTONICITY CHECK ===");
    foreach (var r in allResults.Where(r => r.LensState != null))
    {
        var amps = r.LensState!.Amplitudes;
        bool mono = true;
        for (int i = 1; i < amps.Length; i++)
            if (amps[i] < amps[i - 1] - 0.001) { mono = false; break; }
        scalingWriter.WriteLine($"{r.Label}: {(mono ? "monotonic" : "NOT monotonic")} [{string.Join(",", amps.Select(a => a.ToString("F3")))}]");
    }

    scalingWriter.WriteLine("\n=== ACCESSIBILITY BOUNDARY ===");
    foreach (var r in allResults)
    {
        string acc = r.SEFrobRatios.Length > 1
            ? (r.SEFrobRatios[1] < 0.01 ? "inaccessible" : $"ACCESSIBLE (SE={r.SEFrobRatios[1]:F4})")
            : "n/a";
        scalingWriter.WriteLine($"{r.Label}: 2nd mode {acc} (Re={r.SecondModeRate:F4})");
    }

    // ---- JSON output ----
    using var jsonWriter = new StreamWriter(jsonPath);
    jsonWriter.WriteLine("[");
    for (int i = 0; i < allResults.Count; i++)
    {
        var r = allResults[i];
        jsonWriter.WriteLine("  {");
        jsonWriter.WriteLine($"    \"label\": \"{r.Label}\",");
        jsonWriter.WriteLine($"    \"N\": {r.N},");
        jsonWriter.WriteLine($"    \"topology\": \"{r.Topology}\",");
        jsonWriter.WriteLine($"    \"gammas\": [{string.Join(",", r.Gammas.Select(g => g.ToString("G6", System.Globalization.CultureInfo.InvariantCulture)))}],");
        if (r.SlowModes.Length > 0)
        {
            jsonWriter.WriteLine($"    \"slow_mode_rate\": {(-r.SlowModes[0].Eigenvalue.Real).ToString("G8", System.Globalization.CultureInfo.InvariantCulture)},");
            jsonWriter.WriteLine($"    \"slow_mode_im\": {r.SlowModes[0].Eigenvalue.Imaginary.ToString("G8", System.Globalization.CultureInfo.InvariantCulture)},");
        }
        jsonWriter.WriteLine($"    \"se_frob_ratio\": {(r.SEFrobRatios.Length > 0 ? r.SEFrobRatios[0].ToString("G8", System.Globalization.CultureInfo.InvariantCulture) : "null")},");
        if (r.LensState != null)
        {
            jsonWriter.WriteLine($"    \"psi_opt\": [{string.Join(",", r.LensState.Amplitudes.Select(a => a.ToString("G8", System.Globalization.CultureInfo.InvariantCulture)))}],");
            jsonWriter.WriteLine($"    \"c_slow\": {r.LensState.SlowModeProjection.ToString("G8", System.Globalization.CultureInfo.InvariantCulture)},");
        }
        jsonWriter.WriteLine($"    \"second_mode_rate\": {(double.IsNaN(r.SecondModeRate) ? "null" : r.SecondModeRate.ToString("G8", System.Globalization.CultureInfo.InvariantCulture))},");
        jsonWriter.WriteLine($"    \"second_mode_accessible\": {r.SecondModeAccessible.ToString().ToLower()}");
        jsonWriter.Write("  }");
        if (i < allResults.Count - 1) jsonWriter.Write(",");
        jsonWriter.WriteLine();
    }
    jsonWriter.WriteLine("]");

    Console.WriteLine($"\nTotal sweep time: {sw.Elapsed.TotalMinutes:F1} min");
    Console.WriteLine($"Results in: {outDir}");
    Console.WriteLine("=== LENS SURVEY COMPLETE ===");
}

// ============================================================
// PTF EXPORT: full dense eigendecomp + left & right eigenvectors (N=7)
// ============================================================
// Outputs raw little-endian complex128 binary files for Python.
// Consumer: simulations/eq014_step23_biorth.py
static void RunPtfExport(string resultsDir)
{
    Directory.CreateDirectory(resultsDir);
    const int N = 7;
    const double J = 1.0;
    const double gamma = 0.05;
    int d = 1 << N;
    int d2 = d * d;

    Console.WriteLine("=== PTF DENSE EIGENDECOMP (EQ-014) ===");
    Console.WriteLine($"Date: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
    Console.WriteLine($"Parameters: N={N}, J={J} (XY PTF conv.: H = (J/2)(XX+YY)), gamma={gamma}, Chain, d2={d2}");
    Console.WriteLine($"Per-matrix size: {(long)d2 * d2 * 16 / 1e9:F2} GB (values + 2 eigvec blocks + workspace ~16 GB peak)");
    Console.WriteLine();

    MklDirect.ConfigureThreads(Environment.ProcessorCount, m => Console.WriteLine($"  {m}"));
    Console.WriteLine();

    var sw = new Stopwatch();
    sw.Start();

    // PTF uses XY-only chain with H = Σ (J/2)(X_i X_{i+1} + Y_i Y_{i+1}).
    // Topology.ChainXY halves J and filters PauliTypes to {X, Y} so that
    // BuildHamiltonian yields the correct PTF convention.
    var bonds = Topology.ChainXY(N, Enumerable.Repeat(J, N - 1).ToArray());
    var gammas = Enumerable.Repeat(gamma, N).ToArray();

    Console.WriteLine("[1/4] Building Liouvillian (N=7 uniform XY chain, PTF convention)...");
    sw.Restart();
    var rawL = Liouvillian.BuildDirectRaw(N, bonds, gammas, s => Console.WriteLine($"  {s}"));
    Console.WriteLine($"  Built in {sw.Elapsed.TotalSeconds:F1} s");

    Console.WriteLine("\n[2/4] Dense eigendecomposition (left + right eigenvectors)...");
    sw.Restart();
    var (values, leftVecs, rightVecs) = Liouvillian.GetAllEigenvaluesLeftRightMklRaw(
        rawL, d2, s => Console.WriteLine($"  {s}"));
    Console.WriteLine($"  Eigendecomposition in {sw.Elapsed.TotalMinutes:F1} min");
    rawL = null!;
    GC.Collect();

    // Validation: eigenvalue count, stationary count, palindrome pairing
    Console.WriteLine("\n[3/4] Validation:");
    Console.WriteLine($"  Eigenvalue count: {values.Length} (expected {d2})");
    int stationary = 0;
    double maxAbsIm = 0;
    double mostNegRe = 0;
    foreach (var v in values)
    {
        if (v.Magnitude < 1e-10) stationary++;
        if (Math.Abs(v.Imaginary) > maxAbsIm) maxAbsIm = Math.Abs(v.Imaginary);
        if (v.Real < mostNegRe) mostNegRe = v.Real;
    }
    Console.WriteLine($"  Stationary modes (|λ|<1e-10): {stationary} (F4 expects N+1 = {N + 1})");
    Console.WriteLine($"  Max |Im(λ)|: {maxAbsIm:F6}");
    Console.WriteLine($"  Most negative Re(λ): {mostNegRe:F6} (boundary formula: -2(N-1)γ = {-2.0 * (N - 1) * gamma:F6})");

    // Sum of real parts: palindrome implies sum = -N*γ*d2 (each eigenvalue and its mirror partner sum to -2N*γ)
    double sumRe = 0;
    foreach (var v in values) sumRe += v.Real;
    double expectedSumRe = -(double)d2 * N * gamma;
    Console.WriteLine($"  Σ Re(λ): {sumRe:F3} (palindrome expects {expectedSumRe:F3}, i.e. avg = -N*γ = {-N * gamma:F3})");

    // Biorthogonality spot check: first 5 slow modes
    Console.WriteLine("\n  Biorthogonality spot check (<L_i | R_j> for a few pairs):");
    var idxBySlow = Enumerable.Range(0, values.Length)
        .OrderBy(i => Math.Abs(values[i].Real))
        .Take(5).ToArray();
    for (int a = 0; a < idxBySlow.Length; a++)
    {
        for (int b = 0; b < idxBySlow.Length; b++)
        {
            int i = idxBySlow[a], j = idxBySlow[b];
            Complex dot = Complex.Zero;
            for (int k = 0; k < d2; k++)
            {
                var li = leftVecs[(long)i * d2 + k];
                var rj = rightVecs[(long)j * d2 + k];
                dot += Complex.Conjugate(li) * rj;
            }
            Console.WriteLine($"    <L_{i}|R_{j}> = ({dot.Real:F6},{dot.Imaginary:F6})  (λ_i={values[i].Real:F4}, λ_j={values[j].Real:F4})");
        }
    }

    // Write outputs
    Console.WriteLine("\n[4/4] Writing outputs...");
    sw.Restart();

    string valPath = Path.Combine(resultsDir, "eq014_eigvals_n7.bin");
    string rightPath = Path.Combine(resultsDir, "eq014_right_eigvecs_n7.bin");
    string leftPath = Path.Combine(resultsDir, "eq014_left_eigvecs_n7.bin");
    string metaPath = Path.Combine(resultsDir, "eq014_metadata.json");

    WriteComplexArrayRaw(valPath, values);
    Console.WriteLine($"  Wrote {valPath} ({new FileInfo(valPath).Length / 1024.0:F1} KB)");

    WriteComplexArrayRaw(rightPath, rightVecs);
    Console.WriteLine($"  Wrote {rightPath} ({new FileInfo(rightPath).Length / 1e9:F2} GB)");

    WriteComplexArrayRaw(leftPath, leftVecs);
    Console.WriteLine($"  Wrote {leftPath} ({new FileInfo(leftPath).Length / 1e9:F2} GB)");

    // Metadata JSON
    using (var jw = new StreamWriter(metaPath))
    {
        jw.WriteLine("{");
        jw.WriteLine($"  \"N\": {N},");
        jw.WriteLine($"  \"d\": {d},");
        jw.WriteLine($"  \"d2\": {d2},");
        jw.WriteLine($"  \"J\": {J.ToString("G8", System.Globalization.CultureInfo.InvariantCulture)},");
        jw.WriteLine($"  \"gamma\": {gamma.ToString("G8", System.Globalization.CultureInfo.InvariantCulture)},");
        jw.WriteLine($"  \"topology\": \"chainXY\",");
        jw.WriteLine($"  \"hamiltonian\": \"H = sum_i (J/2)(X_i X_{{i+1}} + Y_i Y_{{i+1}}), PTF convention\",");
        jw.WriteLine($"  \"eigenvalue_count\": {values.Length},");
        jw.WriteLine($"  \"stationary_count\": {stationary},");
        jw.WriteLine($"  \"stationary_expected\": {N + 1},");
        jw.WriteLine($"  \"palindrome_center_re\": {(-N * gamma).ToString("G8", System.Globalization.CultureInfo.InvariantCulture)},");
        jw.WriteLine($"  \"sum_re\": {sumRe.ToString("G8", System.Globalization.CultureInfo.InvariantCulture)},");
        jw.WriteLine($"  \"most_negative_re\": {mostNegRe.ToString("G8", System.Globalization.CultureInfo.InvariantCulture)},");
        jw.WriteLine($"  \"layout\": \"column-major, Fortran-order\",");
        jw.WriteLine($"  \"dtype\": \"complex128 (little-endian)\",");
        jw.WriteLine($"  \"right_eigvecs_convention\": \"L * r_j = lambda_j * r_j ; column j of rightVecs is r_j\",");
        jw.WriteLine($"  \"left_eigvecs_convention\": \"u_j^H * L = lambda_j * u_j^H ; column j of leftVecs is u_j ; biorth dot = conj(u_i)^T * r_j\",");
        jw.WriteLine($"  \"date\": \"{DateTime.Now:yyyy-MM-dd HH:mm:ss}\"");
        jw.WriteLine("}");
    }
    Console.WriteLine($"  Wrote {metaPath}");
    Console.WriteLine($"  Write time: {sw.Elapsed.TotalSeconds:F1} s");

    Console.WriteLine("\n=== PTF EXPORT COMPLETE ===");
}

static void WriteComplexArrayRaw(string path, Complex[] array)
{
    using var fs = new FileStream(path, FileMode.Create, FileAccess.Write, FileShare.None,
        bufferSize: 1 << 20 /*1 MB*/, useAsync: false);
    const int chunkElems = 1 << 20; // 1M Complex = 16 MB per chunk
    var buffer = new byte[chunkElems * 16];
    long total = array.LongLength;
    long written = 0;
    while (written < total)
    {
        long remaining = total - written;
        int thisChunk = (int)Math.Min(chunkElems, remaining);
        // Copy: real as 8 bytes, imag as 8 bytes, for each element
        for (int i = 0; i < thisChunk; i++)
        {
            double re = array[written + i].Real;
            double im = array[written + i].Imaginary;
            Buffer.BlockCopy(BitConverter.GetBytes(re), 0, buffer, i * 16, 8);
            Buffer.BlockCopy(BitConverter.GetBytes(im), 0, buffer, i * 16 + 8, 8);
        }
        fs.Write(buffer, 0, thisChunk * 16);
        written += thisChunk;
    }
}
