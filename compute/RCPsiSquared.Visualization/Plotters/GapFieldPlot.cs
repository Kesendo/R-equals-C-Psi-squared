using System.Linq;
using ScottPlot;

namespace RCPsiSquared.Visualization.Plotters;

/// <summary>The F89 octic branch-locus visual: the octic min-gap intensity over the complex-q plane as
/// a cyberpunk + Matrix heatmap (dark to neon green to cyan-white at the branch points), with the
/// exceptional points (EPs) marked in magenta and the diabolic point in cyan. The data is the
/// flashlight's gap field; this is its colour photograph.</summary>
public static class GapFieldPlot
{
    // a centred glow: translucent halos under a bright core, at the true point positions.
    private static void Glow(Plot plot, (double re, double im)[] pts, Color color, double core)
    {
        double[] xs = pts.Select(p => p.re).ToArray();
        double[] ys = pts.Select(p => p.im).ToArray();
        (double mul, byte a)[] halos = { (5.0, 12), (3.3, 22), (2.1, 42), (1.4, 85) };
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
    }

    public static void Save(
        double[,] intensity, double reLo, double reHi, double imLo, double imHi,
        (double re, double im)[] eps, (double re, double im)[] diabolic,
        string title, string outPath, int width = 1700, int height = 1000)
    {
        var plot = new Plot();

        var hm = plot.Add.Heatmap(intensity);
        hm.Colormap = new CyberpunkMatrixColormap();
        hm.Extent = new CoordinateRect(reLo, reHi, imLo, imHi);

        // Each branch point gets an EXPLICIT centred glow at its true q, so the ring always sits on a green
        // (EP) / cyan (diabolic) dot, independent of how the gap field happened to sample that cusp.
        if (diabolic.Length > 0)
        {
            Glow(plot, diabolic, Color.FromHex("#00fff7"), 7);
            var d = plot.Add.Markers(diabolic.Select(e => e.re).ToArray(), diabolic.Select(e => e.im).ToArray());
            d.MarkerStyle.Shape = MarkerShape.OpenDiamond;
            d.MarkerStyle.Size = 16;
            d.MarkerStyle.LineColor = Color.FromHex("#00fff7");
            d.MarkerStyle.LineWidth = 2;
            d.LegendText = "diabolic point (silent)";
        }

        if (eps.Length > 0)
        {
            Glow(plot, eps, Color.FromHex("#39ff14"), 6);
            var m = plot.Add.Markers(eps.Select(e => e.re).ToArray(), eps.Select(e => e.im).ToArray());
            m.MarkerStyle.Shape = MarkerShape.OpenCircle;
            m.MarkerStyle.Size = 12;
            m.MarkerStyle.LineColor = Color.FromHex("#ff2bd6");
            m.MarkerStyle.LineWidth = 2;
            m.LegendText = "exceptional points (transpositions)";
        }

        plot.FigureBackground.Color = Color.FromHex("#05060a");
        plot.DataBackground.Color = Color.FromHex("#000500");
        plot.Axes.Color(Color.FromHex("#39ff77"));
        plot.Grid.MajorLineColor = Color.FromHex("#0c2a16");
        plot.Title(title);
        plot.XLabel("Re q   (q = J / γ)");
        plot.YLabel("Im q");
        plot.ShowLegend();

        plot.SavePng(outPath, width, height);
    }
}
