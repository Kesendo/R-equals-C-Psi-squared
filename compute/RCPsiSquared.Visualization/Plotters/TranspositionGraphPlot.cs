using System.Collections.Generic;
using System.Linq;
using ScottPlot;

namespace RCPsiSquared.Visualization.Plotters;

/// <summary>The minimal C-T witness for F89_MONODROMY_MIRROR: σ_T is not a symmetry of the braiding.
/// Only the few strands of one witness are drawn (the rest of the connected braid graph, which says
/// monodromy = S_8, is the OTHER result, shown by --root galoismonodromy). σ_T mirrors the strands
/// left-right across the Re λ = −4 fold (twins at ±x, fixed strands on the axis). One real braid edge is
/// drawn solid green; its σ_T-image is drawn dotted, and it is NOT a braid edge. A real connection reflects
/// to a non-connection: σ_T does not commute with the monodromy (it is non-central), so the fold is not a
/// loop-independent symmetry of the braiding.</summary>
public static class TranspositionGraphPlot
{
    /// <param name="xs">node x (−1 left twin, 0 on the fold, +1 right twin)</param>
    /// <param name="ys">node y (schematic rows)</param>
    /// <param name="labels">strand index per node</param>
    /// <param name="onFold">true if the strand sits on the Re λ = −4 fold (σ_T-fixed)</param>
    /// <param name="braid">the real braid edge, as node indices into xs/ys</param>
    /// <param name="ghost">its σ_T-image, as node indices (NOT a braid edge)</param>
    public static void SaveWitness(
        double[] xs, double[] ys, int[] labels, bool[] onFold,
        (int a, int b) braid, (int a, int b) ghost,
        string title, string outPath, int width = 1200, int height = 950)
    {
        var plot = new Plot();

        // the σ_T fold axis (the Re λ = −4 mirror): vertical dashed amber.
        var ax = plot.Add.VerticalLine(0);
        ax.Color = Color.FromHex("#ffe700");
        ax.LinePattern = LinePattern.Dashed;
        ax.LineWidth = 1.5f;
        ax.LegendText = "the Re λ = −4 fold  (σ_T = left-right mirror)";

        // the real braid edge: solid bright green ("this IS a braid").
        var be = plot.Add.Scatter(new[] { xs[braid.a], xs[braid.b] }, new[] { ys[braid.a], ys[braid.b] });
        be.Color = Color.FromHex("#39ff14");
        be.LineWidth = 6;
        be.MarkerSize = 0;
        be.LegendText = $"a braid edge ({labels[braid.a]} {labels[braid.b]})  —  an EP transposition";

        // its σ_T-image: dotted magenta ("σ_T reflects it here, but this is NOT a braid").
        var ge = plot.Add.Scatter(new[] { xs[ghost.a], xs[ghost.b] }, new[] { ys[ghost.a], ys[ghost.b] });
        ge.Color = Color.FromHex("#ff2bd6");
        ge.LinePattern = LinePattern.Dotted;
        ge.LineWidth = 5;
        ge.MarkerSize = 0;
        ge.LegendText = $"its σ_T-image ({labels[ghost.a]} {labels[ghost.b]})  —  NOT a braid edge";

        // nodes: gold on the fold (Re λ = −4), cyan off it (the mirror-twins).
        Glow(plot, Pts(xs, ys, onFold, true), Color.FromHex("#ffae00"), 20, "a strand on the fold (Re λ = −4)");
        Glow(plot, Pts(xs, ys, onFold, false), Color.FromHex("#08f7fe"), 20, "mirror-twin strands (σ_T swaps them)");

        for (int i = 0; i < xs.Length; i++)
        {
            var t = plot.Add.Text(labels[i].ToString(), xs[i], ys[i]);
            t.LabelFontColor = Color.FromHex("#05060a");
            t.LabelFontSize = 18;
            t.LabelBold = true;
            t.Alignment = Alignment.MiddleCenter;
        }

        double yMin = ys.Min(), yMax = ys.Max();
        plot.Axes.SetLimits(-2.0, 2.0, yMin - 0.8, yMax + 0.8);
        plot.Axes.Frameless();                                  // schematic: no quantitative axes
        plot.FigureBackground.Color = Color.FromHex("#05060a");
        plot.DataBackground.Color = Color.FromHex("#05060a");
        plot.Title(title);
        plot.ShowLegend(Edge.Bottom);

        plot.SavePng(outPath, width, height);
    }

    private static (double re, double im)[] Pts(double[] xs, double[] ys, bool[] onFold, bool wantFold)
        => Enumerable.Range(0, xs.Length).Where(i => onFold[i] == wantFold).Select(i => (xs[i], ys[i])).ToArray();

    private static void Glow(Plot plot, (double re, double im)[] pts, Color color, double core, string legend)
    {
        if (pts.Length == 0) return;
        double[] xs = pts.Select(p => p.re).ToArray();
        double[] ys = pts.Select(p => p.im).ToArray();
        (double mul, byte a)[] halos = { (4.0, 18), (2.7, 32), (1.8, 65) };
        foreach (var (mul, a) in halos)
        {
            var h = plot.Add.Markers(xs, ys);
            h.MarkerStyle.Shape = MarkerShape.FilledCircle;
            h.MarkerStyle.Size = (float)(core * mul);
            h.MarkerStyle.FillColor = color.WithAlpha(a);
            h.MarkerStyle.LineWidth = 0;
        }
        var c = plot.Add.Markers(xs, ys);
        c.MarkerStyle.Shape = MarkerShape.FilledCircle;
        c.MarkerStyle.Size = (float)core;
        c.MarkerStyle.FillColor = color;
        c.MarkerStyle.LineWidth = 0;
        c.LegendText = legend;
    }
}
