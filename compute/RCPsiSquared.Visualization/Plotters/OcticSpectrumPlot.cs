using System.Linq;
using ScottPlot;

namespace RCPsiSquared.Visualization.Plotters;

/// <summary>The F89 path-3 octic spectrum in the λ-plane, time-killed (each mode = one point λ), drawn
/// in the symphony reel palette so it lays directly beside the molecule-spectrum images of
/// experiments/THE_SHARED_SKELETON.md. The AT-locked roots (cyan) sit on the absorption rungs
/// Re λ = −2γ, −6γ; the H_B-mixed octic (pink, Galois S_8) spreads between them. Re λ is the dephasing
/// diagonal (the watched axis z), Im λ the Hamiltonian's motion (x, y); the palindrome centre is Re = −4.</summary>
public static class OcticSpectrumPlot
{
    private static readonly Color Bg = Color.FromHex("#080b12");
    private static readonly Color Panel = Color.FromHex("#0c1018");
    private static readonly Color Fg = Color.FromHex("#7ef9ff");
    private static readonly Color GridCol = Color.FromHex("#15303a");
    private static readonly Color NeonGreen = Color.FromHex("#39ff14");
    private static readonly Color Pink = Color.FromHex("#fe53bb");
    private static readonly Color Cyan = Color.FromHex("#08f7fe");
    private static readonly Color Amber = Color.FromHex("#ffe700");

    // a glow scatter: several translucent halos under a bright core, mimicking the reel's glow_scatter.
    private static void Glow(Plot plot, (double re, double im)[] pts, Color color, double core, string legend)
    {
        double[] xs = pts.Select(p => p.re).ToArray();
        double[] ys = pts.Select(p => p.im).ToArray();
        (double mul, byte a)[] halos = { (5.0, 14), (3.4, 24), (2.2, 45), (1.45, 90) };
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

    public static void Save((double re, double im)[] octic, (double re, double im)[] at,
        double centerRe, double q, string title, string outPath, int width = 1400, int height = 1000)
    {
        var plot = new Plot();

        var cl = plot.Add.VerticalLine(centerRe);
        cl.Color = Amber;
        cl.LinePattern = LinePattern.Dashed;
        cl.LineWidth = 1;
        cl.LegendText = $"palindrome centre  Re λ = {centerRe:0.#}";

        Glow(plot, at, Cyan, 11, "AT-locked (free-fermion; rungs −2γ, −6γ)");
        Glow(plot, octic, Pink, 11, "H_B-mixed octic (Galois S_8, no radical closure)");

        plot.FigureBackground.Color = Bg;
        plot.DataBackground.Color = Panel;
        plot.Axes.Color(Fg);
        plot.Grid.MajorLineColor = GridCol;
        plot.Title(title);
        plot.XLabel("Re λ / γ   (dephasing diagonal, the watched axis z)");
        plot.YLabel("Im λ   (Hamiltonian motion, the plane x, y)");
        plot.ShowLegend();

        plot.SavePng(outPath, width, height);
    }
}
