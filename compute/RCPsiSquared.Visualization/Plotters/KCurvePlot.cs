using RCPsiSquared.Core.Resonance;
using ScottPlot;

namespace RCPsiSquared.Visualization.Plotters;

/// <summary>Plot raw K_b(Q) curves per bond, plus the bond-class averages.</summary>
public static class KCurvePlot
{
    public static Plot Build(KCurve curve, PlotConfig? config = null)
    {
        var plot = new Plot();

        double[] qGrid = curve.QGrid.ToArray();

        // Per-bond curves in light grey
        for (int b = 0; b < curve.NumBonds; b++)
        {
            var per = plot.Add.Scatter(qGrid, curve.BondCurve(b));
            per.Color = ScottPlot.Colors.LightGray;
            per.LineWidth = 1f;
            per.MarkerSize = 0;
            per.LegendText = $"bond {b}";
        }

        // Bond-class averages on top
        var interior = plot.Add.Scatter(qGrid, curve.BondClassAverage(BondClass.Interior));
        interior.Color = BondClassColors.Interior;
        interior.LineWidth = 2f;
        interior.MarkerSize = 0;
        interior.LegendText = "Interior (mean)";

        var endpoint = plot.Add.Scatter(qGrid, curve.BondClassAverage(BondClass.Endpoint));
        endpoint.Color = BondClassColors.Endpoint;
        endpoint.LineWidth = 2f;
        endpoint.MarkerSize = 0;
        endpoint.LegendText = "Endpoint (mean)";

        plot.Title($"K_CC_pr(Q) — N={curve.Block.N}, n={curve.Block.LowerPopcount}, γ₀={curve.Block.GammaZero}");
        plot.XLabel("Q = J/γ₀");
        plot.YLabel("|K|max over t");
        plot.ShowLegend();

        return plot;
    }

    public static void Save(KCurve curve, string outputPath, PlotConfig? config = null)
    {
        var c = config ?? PlotConfig.Default;
        var plot = Build(curve, c);
        plot.SavePng(outputPath, c.Width, c.Height);
    }
}
