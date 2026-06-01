using System.Globalization;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Cli.Commands;

/// <summary>The assembly read across a Q-sweep, the whole locked picture played through on one
/// substrate: for each Q the slowest mode's rate, its drain depth (= light n_XY = rate), its parity
/// rung (flow / birth), the Absorption-Theorem cross-check (rate = 2·Σγ·light), and the per-site
/// carrier vector; then once the birth-canal vs sterile verdict and the maximal-saturation ceiling.
/// Chain (open line) or ring (the aromatic substrate: benzene C₆, cyclobutadiene C₄).</summary>
public static class AssemblyCommand
{
    public static int Run(string[] args)
    {
        var p = new ArgParser(args);
        p.RequireNoPositional();
        int N = p.RequireInt("N");

        string qListStr = p.OptionalString("q-list") ?? "0.5,1.0,1.5,2.5,5.0,20.0";
        var qGrid = qListStr.Split(',')
            .Select(s => double.Parse(s.Trim(), CultureInfo.InvariantCulture)).ToArray();

        FlowTopology topology = (p.OptionalString("topology") ?? "chain").ToLowerInvariant() switch
        {
            "chain" => FlowTopology.Chain,
            "ring" => FlowTopology.Ring,
            var other => throw new ArgumentException($"unknown topology: {other} (chain|ring)"),
        };

        double[]? profile = null;
        string? profileStr = p.OptionalString("gamma-profile");
        if (profileStr is not null)
        {
            profile = profileStr.Split(',')
                .Select(s => double.Parse(s.Trim(), CultureInfo.InvariantCulture)).ToArray();
            if (p.HasFlag("fix-total")) profile = PostEpFlowField.NormalizeToTotal(profile, N);
        }

        // A coarse τ-grid is enough for the between-block saturation reading.
        var tau = new double[40];
        for (int i = 0; i < tau.Length; i++) tau[i] = 6.0 * i / (tau.Length - 1);
        var field = new PostEpFlowField(N, qGrid, tau, profile, topology);

        string profLabel = profile is null
            ? "uniform"
            : "[" + string.Join(",", profile.Select(w => w.ToString("0.###", CultureInfo.InvariantCulture))) + "]";
        Console.WriteLine($"Assembly read: N={N}, topology={topology}, gamma-profile={profLabel}");
        Console.WriteLine($"  {"Q",7}  {"rate",8}  {"depth",6}  {"rung",11}  {"absorp",7}  per-site light <X/Y at k>");
        foreach (double q in qGrid)
        {
            var a = field.ReadAssembly(q);
            string rung = a.OnBirthRail ? "BIRTH(odd)" : "flow(even)";
            string ok = Math.Abs(a.SlowestRate - a.AbsorptionRate) < 1e-6 ? "ok" : "MISMATCH";
            string light = "[" + string.Join(" ",
                a.PerSiteLight.Select(x => x.ToString("0.00", CultureInfo.InvariantCulture))) + "]";
            Console.WriteLine($"  {q,7:0.##}  {a.SlowestRate,8:0.0000}  {a.SlowestDepth,6:0.000}  {rung,11}  {ok,7}  deg{a.Degeneracy,3}  {light}");
        }

        Console.WriteLine();
        string zone = field.IsInBirthCanal
            ? $"BIRTH CANAL (Q-drift δ = {field.BirthCanalDeviation.ToString("0.0000", CultureInfo.InvariantCulture)})"
            : "STERILE (the birth channel's lifetime is Q-independent)";
        Console.WriteLine($"  zone: {zone}");
        Console.WriteLine($"  max saturation ceiling C_block = {PostEpFlowField.MaxSaturationCeiling.ToString("0.00", CultureInfo.InvariantCulture)} (= 1/4); " +
            $"flow's peak between-block = {field.PeakBetweenBlockSaturation.ToString("E2", CultureInfo.InvariantCulture)} (~0, the flow rides the even rail)");
        return 0;
    }
}
