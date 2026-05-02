using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Resonance;
using RCPsiSquared.Visualization.Plotters;

namespace RCPsiSquared.Cli.Commands;

public enum PlotKind
{
    KCurve,
    Shape,
}

public static class PlotCommand
{
    public static int Run(string[] args)
    {
        if (args.Length == 0)
        {
            Console.Error.WriteLine($"plot needs a kind: {string.Join(" | ", Enum.GetNames<PlotKind>()).ToLowerInvariant()}");
            return 2;
        }
        if (!Enum.TryParse<PlotKind>(args[0], ignoreCase: true, out var kind))
        {
            Console.Error.WriteLine($"unknown plot kind: {args[0]}");
            return 2;
        }

        var p = new ArgParser(args[1..]);
        p.RequireNoPositional();
        int N = p.RequireInt("N");
        int n = p.RequireInt("n");
        double gamma = p.RequireDouble("gamma");
        string outPath = p.OptionalString("out") ?? "plot.png";

        var block = new CoherenceBlock(N, n, gamma);
        var scan = new ResonanceScan(block);
        var curve = scan.ComputeKCurve();

        switch (kind)
        {
            case PlotKind.KCurve:
                KCurvePlot.Save(curve, outPath);
                break;
            case PlotKind.Shape:
                ResonanceShapePlot.Save(new[] { ($"c{block.C}_N{N}", curve) }, outPath);
                break;
        }

        Console.WriteLine($"saved {outPath}");
        return 0;
    }
}
