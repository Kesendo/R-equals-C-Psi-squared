using System.Text.Json;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Cli.Commands;

public static class ScanCommand
{
    public static int Run(string[] args)
    {
        var p = new ArgParser(args);
        p.RequireNoPositional();
        int N = p.RequireInt("N");
        int n = p.RequireInt("n");
        double gamma = p.RequireDouble("gamma");
        string? outPath = p.OptionalString("out");

        var block = new CoherenceBlock(N, n, gamma);
        var scan = new ResonanceScan(block);
        var curve = scan.ComputeKCurve();

        var interior = curve.Peak(BondClass.Interior);
        var endpoint = curve.Peak(BondClass.Endpoint);

        var result = new
        {
            input = new { N, lower_popcount = n, gamma_0 = gamma },
            block_dim = block.Basis.MTotal,
            chromaticity = block.C,
            num_bonds = block.NumBonds,
            interior = SummarisePeak(interior),
            endpoint = SummarisePeak(endpoint),
        };

        var json = JsonSerializer.Serialize(result, new JsonSerializerOptions { WriteIndented = true });
        if (outPath is null) Console.WriteLine(json);
        else File.WriteAllText(outPath, json);
        return 0;
    }

    private static object SummarisePeak(PeakResult peak) => new
    {
        q_peak = peak.QPeak,
        k_max = peak.KMax,
        hwhm_left = peak.HwhmLeft,
        hwhm_right = peak.HwhmRight,
        hwhm_left_over_q_peak = peak.HwhmLeftOverQPeak,
        asymmetry = peak.Asymmetry,
    };
}
