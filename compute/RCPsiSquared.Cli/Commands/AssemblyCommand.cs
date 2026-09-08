using System.Globalization;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Ptf;

namespace RCPsiSquared.Cli.Commands;

/// <summary>The assembly read across a Q-sweep, the whole locked picture played through on one
/// substrate: for each canonical Q the global slowest mode's rate, its dimensionless drain depth
/// (= light n_XY), its projector parity support, the distinct weighted Absorption rate
/// 2·Σγ_l·light_l, and the per-site carrier vector; then, where scoped, the N=5 open-chain
/// birth-canal vs rate-sterile verdict and the maximal-saturation ceiling.
/// Chain (open line) or ring (the aromatic substrate: benzene C₆, cyclobutadiene C₄).</summary>
public static class AssemblyCommand
{
    public static int Run(string[] args)
    {
        var p = new ArgParser(args);
        p.RequireNoPositional();
        int N = p.RequireInt("N");

        string qListStr = p.OptionalString("q-list") ?? "1.0,2.0,3.0,5.0,10.0,40.0";
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
        Console.WriteLine($"  {"Q",7}  {"edge",8}  {"mean",8}  {"depth",6}  {"parity",13}  {"abs(mean)",9}  per-site light <X/Y at k>");
        foreach (double q in qGrid)
        {
            var a = field.ReadAssembly(q);
            string rung = a.Parity switch
            {
                SlowLightParity.Odd => "odd support",
                SlowLightParity.Even => "even support",
                _ => "mixed support",
            };
            string ok = a.AbsorptionResidual < 1e-6 ? "ok" : "MISMATCH";
            string light = "[" + string.Join(" ",
                a.PerSiteLight.Select(x => x.ToString("0.00", CultureInfo.InvariantCulture))) + "]";
            Console.WriteLine($"  {q,7:0.##}  {a.SlowestRate,8:0.0000}  {a.SlowClusterMeanRate,8:0.0000}  {a.SlowestDepth,6:0.000}  {rung,13}  {ok,9}  cluster{a.SlowClusterDimension,3}  {light}");
        }

        Console.WriteLine();
        string zone = field.HasBirthCanalClassification
            ? field.IsInBirthCanal
                ? $"BIRTH CANAL (Q-drift δ = {field.BirthCanalDeviation.ToString("0.0000", CultureInfo.InvariantCulture)})"
                : "RATE-STERILE on the verified N=5 surface"
            : $"UNCLASSIFIED outside the N=5 surface (global slow-rate drift δ = " +
              $"{field.GlobalSlowestRateDeviation.ToString("0.0000", CultureInfo.InvariantCulture)})";
        Console.WriteLine($"  zone: {zone}");
        Console.WriteLine($"  max saturation ceiling C_block = {PostEpFlowField.MaxSaturationCeiling.ToString("0.00", CultureInfo.InvariantCulture)} (= 1/4); " +
            $"flow's peak between-block = {field.PeakBetweenBlockSaturation.ToString("E2", CultureInfo.InvariantCulture)} " +
            "(~0, the prepared single-excitation flow has fixed excitation number)");
        return 0;
    }
}
