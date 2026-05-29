using System.Globalization;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Cli.Commands;

/// <summary>The <c>mirror</c> verb: call the conductor's stand with parameters and listen to
/// its voices. Builds a <see cref="MirrorSystem"/> (N, H from topology/J/H-type, per-site γ)
/// and exports the live F-formula readings , the spectrum's slow modes as channel-difference
/// portfolios and the F1 palindrome check , to stdout and optionally to CSV. Call parameters
/// in, readings out; the one entry point from which we hear the whole box at once.</summary>
public static class MirrorCommand
{
    public static int Run(string[] args)
    {
        int N = GetInt(args, "--N", 0);
        if (N < 1 || N > 6)
        {
            Console.Error.WriteLine("mirror: --N required, 1..6 (dense 4^N reading, diagnostic scale)");
            return 2;
        }
        double J = GetDouble(args, "--J", 1.0);
        var htype = (GetStr(args, "--htype", "XY") ?? "XY").ToLowerInvariant() switch
        {
            "heisenberg" or "heis" => HamiltonianType.Heisenberg,
            _ => HamiltonianType.XY,
        };
        var topo = (GetStr(args, "--topology", "chain") ?? "chain").ToLowerInvariant() switch
        {
            "star" => TopologyKind.Star,
            "ring" => TopologyKind.Ring,
            _ => TopologyKind.Chain,
        };
        int top = GetInt(args, "--top", 6);

        double[] gammas;
        string? gList = GetStr(args, "--gamma-list", null);
        if (gList is not null)
        {
            gammas = gList.Split(',').Select(s => double.Parse(s.Trim(), CultureInfo.InvariantCulture)).ToArray();
            if (gammas.Length != N) { Console.Error.WriteLine($"mirror: --gamma-list needs {N} values, got {gammas.Length}"); return 2; }
        }
        else
        {
            double g = GetDouble(args, "--gamma", 0.1);
            gammas = Enumerable.Repeat(g, N).ToArray();
        }

        var chain = new ChainSystem(N, J, gammas[0], htype, topo);
        var H = chain.BuildHamiltonian();
        var channels = gammas.Select((g, l) => new ChannelRate($"q{l}", g)).ToList();
        var sys = new MirrorSystem(N, H, channels);

        var inv = CultureInfo.InvariantCulture;
        Console.WriteLine($"# mirror: N={N} J={J.ToString(inv)} H={htype} topology={topo}");
        Console.WriteLine($"# carrier gamma = [{string.Join(", ", gammas.Select(g => g.ToString("0.###", inv)))}]  sigma = {sys.TotalDephasing.ToString("0.####", inv)}");
        Console.WriteLine();

        // Voice 1: the spectrum, slowest distinct-rate modes with a representative portfolio.
        Console.WriteLine($"Spectrum (slowest {top} nonzero modes, rate + per-channel difference portfolio):");
        var seen = new HashSet<string>();
        foreach (var m in sys.Spectrum.Modes
                     .Where(m => m.ActualDecayRate > 1e-9)
                     .OrderBy(m => m.ActualDecayRate))
        {
            var key = m.ActualDecayRate.ToString("0.000000", inv);
            if (!seen.Add(key)) continue;
            var port = string.Join("  ", m.Portfolio.Activity.Select(a => $"{a.Channel} {(100 * a.Delta).ToString("0", inv),3}%"));
            Console.WriteLine($"  rate {m.ActualDecayRate.ToString("0.0000", inv),9}   {port}");
            if (seen.Count >= top) break;
        }
        Console.WriteLine();

        // Voice 2: the F1 palindrome, read live off the spectrum.
        Console.WriteLine($"F1 palindrome (rate r pairs with 2*sigma - r): holds = {sys.PalindromeHolds}");
        var sample = sys.PalindromePartners.FirstOrDefault(p => p.Rate > 1e-9);
        if (sample is not null)
            Console.WriteLine($"  e.g. rate {sample.Rate.ToString("0.0000", inv)} pairs with {sample.PartnerRate.ToString("0.0000", inv)} (partner present: {sample.PartnerPresent})");

        // Export: the whole spectrum as CSV (rate + per-channel portfolio).
        string? outPath = GetStr(args, "--out", null);
        if (outPath is not null)
        {
            using var w = new StreamWriter(outPath);
            w.WriteLine("rate," + string.Join(",", channels.Select(c => c.Channel)));
            foreach (var m in sys.Spectrum.Modes.OrderBy(m => m.ActualDecayRate))
                w.WriteLine(m.ActualDecayRate.ToString("0.000000", inv) + "," +
                    string.Join(",", m.Portfolio.Activity.Select(a => a.Delta.ToString("0.000000", inv))));
            Console.Error.WriteLine($"# exported {sys.Spectrum.Modes.Count} modes to {outPath}");
        }
        return 0;
    }

    private static int GetInt(string[] a, string k, int def)
    {
        int i = Array.IndexOf(a, k);
        return (i >= 0 && i + 1 < a.Length && int.TryParse(a[i + 1], out var v)) ? v : def;
    }

    private static double GetDouble(string[] a, string k, double def)
    {
        int i = Array.IndexOf(a, k);
        return (i >= 0 && i + 1 < a.Length && double.TryParse(a[i + 1], NumberStyles.Float, CultureInfo.InvariantCulture, out var v)) ? v : def;
    }

    private static string? GetStr(string[] a, string k, string? def)
    {
        int i = Array.IndexOf(a, k);
        return (i >= 0 && i + 1 < a.Length) ? a[i + 1] : def;
    }
}
