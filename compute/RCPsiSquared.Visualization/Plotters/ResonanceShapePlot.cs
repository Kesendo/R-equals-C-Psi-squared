using RCPsiSquared.Core.Resonance;
using ScottPlot;

namespace RCPsiSquared.Visualization.Plotters;

/// <summary>Universal-shape collapse plot: y = K(Q)/|K|max as function of x = (Q − Q_peak)/Q_peak.
///
/// <para>Visualises Statement 2 of F86: shape universality within bond classes. Multiple
/// (c, N) cases overlay onto two universal shapes (one per class).</para>
/// </summary>
public static class ResonanceShapePlot
{
    public static Plot Build(IReadOnlyList<(string Label, KCurve Curve)> cases, PlotConfig? config = null)
    {
        var plot = new Plot();

        foreach (var (label, curve) in cases)
        {
            AddNormalizedCurve(plot, curve, BondClass.Interior, label, BondClassColors.Interior);
            AddNormalizedCurve(plot, curve, BondClass.Endpoint, label, BondClassColors.Endpoint);
        }

        plot.Title("F86 universal-shape collapse: y = K(Q)/|K|max vs x = (Q−Q*)/Q*");
        plot.XLabel("(Q − Q_peak) / Q_peak");
        plot.YLabel("K(Q) / |K|max");
        plot.ShowLegend();

        return plot;
    }

    private static void AddNormalizedCurve(Plot plot, KCurve curve, BondClass cls, string label, ScottPlot.Color color)
    {
        double[] avg = curve.BondClassAverage(cls);
        var peak = curve.Peak(avg);
        double qStar = peak.QPeak;
        double kMax = peak.KMax;
        if (kMax <= 0) return;

        double[] xs = new double[curve.QGrid.Count];
        double[] ys = new double[curve.QGrid.Count];
        for (int i = 0; i < curve.QGrid.Count; i++)
        {
            xs[i] = (curve.QGrid[i] - qStar) / qStar;
            ys[i] = avg[i] / kMax;
        }
        var s = plot.Add.Scatter(xs, ys);
        s.Color = color;
        s.LineWidth = 1.5f;
        s.MarkerSize = 0;
        s.LegendText = $"{label} {cls}";
    }

    public static void Save(IReadOnlyList<(string Label, KCurve Curve)> cases, string outputPath, PlotConfig? config = null)
    {
        var c = config ?? PlotConfig.Default;
        var plot = Build(cases, c);
        plot.SavePng(outputPath, c.Width, c.Height);
    }
}
