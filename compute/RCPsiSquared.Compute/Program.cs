using System.Diagnostics;
using System.Numerics;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Compute;

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

var outputPath = Path.Combine(
    @"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared",
    "simulations", "results_csharp_compute.txt");

using var writer = new StreamWriter(outputPath);

void Log(string msg)
{
    Console.WriteLine(msg);
    writer.WriteLine(msg);
    writer.Flush();
}

const double gamma = 0.05;
var sw = new Stopwatch();

Log("=" + new string('=', 79));
Log("R=CPsi^2 COMPUTE ENGINE (C#/.NET + MKL)");
Log($"Started: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
Log($"Machine: {Environment.ProcessorCount} cores, {GC.GetGCMemoryInfo().TotalAvailableMemoryBytes / 1e9:F1} GB RAM");
Log("=" + new string('=', 79));

// ============================================================
// BENCHMARK N=2 to N=7
// ============================================================
Log("\n### BENCHMARK: Liouvillian eigendecomposition timing");
Log($"{"N",4} | {"Matrix",8} | {"Build(ms)",10} | {"Eigen(ms)",10} | {"Total(ms)",10} | {"Rates",7} | {"Mirror",7}");
Log(new string('-', 70));

var rateCounts = new Dictionary<int, int>();

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
Console.WriteLine($"\n>>> Results saved to: {outputPath}");
