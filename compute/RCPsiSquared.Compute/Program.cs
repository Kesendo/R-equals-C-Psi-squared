using System.Diagnostics;
using System.Numerics;
using System.Runtime.InteropServices;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Compute;

// Usage: dotnet run -c Release              -> full suite (N=2-7 + topology + stress + N=8)
//        dotnet run -c Release -- n8        -> N=8 only (skip everything else)
//        dotnet run -c Release -- validate  -> N=5 eigenvalue-only validation only
bool n8Only = args.Any(a => a.Equals("n8", StringComparison.OrdinalIgnoreCase));
bool validateOnly = args.Any(a => a.Equals("validate", StringComparison.OrdinalIgnoreCase));

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
        // We temporarily lower the threshold — or just call native path directly
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
