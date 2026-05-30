using System.Globalization;
using System.Numerics;
using System.Text;
using RCPsiSquared.Core.Inspection;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Visualization.Inspection;

/// <summary>Draws an <see cref="InspectablePayload"/> as ASCII art lines, the picture the
/// <see cref="ConsoleTreeRenderer"/> hangs under a node when payload drawing is on. This is the
/// step from naming a leaf's shape (<c>vector[12]</c>) to showing its data:
///
/// <list type="bullet">
///   <item><see cref="InspectablePayload.Vector"/> , one magnitude bar per component (|z|).</item>
///   <item><see cref="InspectablePayload.MatrixView"/> , a magnitude heatmap, the ramp
///         <c> ·░▒▓█</c> by |z|/max.</item>
///   <item><see cref="InspectablePayload.Curve"/> , a multi-row vertical bar plot with a y-axis
///         scale and the x range.</item>
/// </list>
///
/// <para><see cref="InspectablePayload.Real"/> / <see cref="InspectablePayload.Scalar"/> /
/// <see cref="InspectablePayload.None"/> return no lines: a single value reads better inline.
/// Every renderer is width-budgeted and downsamples (a 400-point curve, a 64×64 matrix) so the
/// drawing never exceeds the terminal.</para></summary>
public static class AsciiPayloadPlot
{
    // Output is culture-invariant (dot decimals), like the rest of the CLI / JSON / CSV surface.
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    // Horizontal sub-cell fill for the magnitude bars (1/8 .. 7/8 of a cell); full cell = '█'.
    private static readonly char[] BarEighths = { ' ', '▏', '▎', '▍', '▌', '▋', '▊', '▉' };
    // Heatmap ramp, lightest to densest, indexed by round(|z|/max · 5).
    private static readonly char[] HeatRamp = { ' ', '·', '░', '▒', '▓', '█' };
    // Vertical sub-cell fill for the curve bars (1/8 .. 7/8, bottom-anchored); full cell = '█'.
    private const string ColEighths = "▁▂▃▄▅▆▇";

    /// <summary>Render whichever payload variant carries drawable data; empty for scalar/none.</summary>
    public static IReadOnlyList<string> Render(InspectablePayload payload, int width = 64) => payload switch
    {
        InspectablePayload.Vector v => Vector(v.Values, v.ComponentLabels, width),
        InspectablePayload.MatrixView m => Matrix(m.Values, width),
        InspectablePayload.Curve c => Curve(c.X, c.Y, c.XLabel, width),
        _ => Array.Empty<string>(),
    };

    /// <summary>One row per component: <c>label  bar  |z|</c>, bars scaled to the largest
    /// magnitude. A vector longer than <paramref name="maxRows"/> is reduced to its top rows by
    /// magnitude (in original index order) with a note, so a 4^N eigenvector still fits.</summary>
    public static IReadOnlyList<string> Vector(ComplexVector values, IReadOnlyList<string>? labels = null,
        int width = 64, int maxRows = 16)
    {
        int n = values.Count;
        if (n == 0) return Array.Empty<string>();

        var mag = new double[n];
        for (int i = 0; i < n; i++) mag[i] = values[i].Magnitude;

        bool truncated = n > maxRows;
        int[] shown = truncated
            ? Enumerable.Range(0, n).OrderByDescending(i => mag[i]).Take(maxRows).OrderBy(i => i).ToArray()
            : Enumerable.Range(0, n).ToArray();

        double max = mag.Max();
        int labelW = shown.Select(i => Label(labels, i).Length).DefaultIfEmpty(1).Max();
        int barW = Math.Max(8, width - labelW - 9);

        var lines = new List<string>(shown.Length + 1);
        foreach (int i in shown)
            lines.Add($"{Label(labels, i).PadRight(labelW)} {Bar(mag[i], max, barW)} {mag[i].ToString("0.000", Inv).PadLeft(6)}");
        if (truncated) lines.Add($"(top {maxRows} of {n} by |·|)");
        return lines;
    }

    /// <summary>A magnitude heatmap, one glyph per cell from the ramp <c> ·░▒▓█</c>; a matrix
    /// larger than <paramref name="maxDim"/> (or the width budget) is block-averaged down. The
    /// footer carries the magnitude scale and any downsampling.</summary>
    public static IReadOnlyList<string> Matrix(ComplexMatrix m, int width = 64, int maxDim = 24)
    {
        int R = m.RowCount, C = m.ColumnCount;
        if (R == 0 || C == 0) return Array.Empty<string>();

        int rOut = Math.Min(R, maxDim);
        int cOut = Math.Min(C, Math.Min(maxDim, Math.Max(1, width / 2)));
        var cell = new double[rOut, cOut];
        double max = 0;
        for (int i = 0; i < rOut; i++)
            for (int j = 0; j < cOut; j++)
            {
                int r0 = i * R / rOut, r1 = Math.Max(r0 + 1, (i + 1) * R / rOut);
                int c0 = j * C / cOut, c1 = Math.Max(c0 + 1, (j + 1) * C / cOut);
                double s = 0; int cnt = 0;
                for (int r = r0; r < r1; r++)
                    for (int c = c0; c < c1; c++) { s += m[r, c].Magnitude; cnt++; }
                double v = cnt > 0 ? s / cnt : 0;
                cell[i, j] = v;
                if (v > max) max = v;
            }

        var lines = new List<string>(rOut + 1);
        for (int i = 0; i < rOut; i++)
        {
            var sb = new StringBuilder();
            for (int j = 0; j < cOut; j++) { sb.Append(Glyph(cell[i, j], max)); sb.Append(' '); }
            lines.Add(sb.ToString().TrimEnd());
        }
        string foot = $"|·| max {max.ToString("0.###", Inv)}";
        if (R > rOut || C > cOut) foot += $", {R}×{C}→{rOut}×{cOut}";
        lines.Add(foot);
        return lines;
    }

    /// <summary>A multi-row vertical bar plot: the curve downsampled to one column per cell,
    /// each column a bar whose height is (y−ymin)/(ymax−ymin), with the y-axis min/max labelled
    /// and the x range on the base line.</summary>
    public static IReadOnlyList<string> Curve(IReadOnlyList<double> x, IReadOnlyList<double> y,
        string? xLabel = null, int width = 48, int height = 5)
    {
        int n = y.Count;
        if (n == 0) return Array.Empty<string>();

        int cols = Math.Min(n, width);
        var col = new double[cols];
        for (int j = 0; j < cols; j++)
        {
            int a = j * n / cols, b = Math.Max(a + 1, (j + 1) * n / cols);
            double s = 0; int c = 0;
            for (int k = a; k < b; k++) { s += y[k]; c++; }
            col[j] = c > 0 ? s / c : 0;
        }

        double ymin = col.Min(), ymax = col.Max(), span = ymax - ymin;
        string top = ymax.ToString("0.###", Inv), bot = ymin.ToString("0.###", Inv);
        int axW = Math.Max(top.Length, bot.Length);

        var lines = new List<string>(height + 1);
        for (int row = height; row >= 1; row--)
        {
            var sb = new StringBuilder();
            string lab = row == height ? top.PadLeft(axW) : row == 1 ? bot.PadLeft(axW) : new string(' ', axW);
            sb.Append(lab);
            sb.Append(row == height || row == 1 ? " ┤" : " │");
            foreach (double v in col)
            {
                double frac = span > 1e-300 ? (v - ymin) / span : 0.0;
                sb.Append(Cell(frac * height - (row - 1)));
            }
            lines.Add(sb.ToString());
        }
        string x0 = x.Count > 0 ? x[0].ToString("0.###", Inv) : "";
        string x1 = x.Count > 0 ? x[^1].ToString("0.###", Inv) : "";
        lines.Add($"{new string(' ', axW)} └ {xLabel ?? "x"} {x0} .. {x1}");
        return lines;
    }

    private static string Label(IReadOnlyList<string>? labels, int i) =>
        labels is not null && i < labels.Count ? labels[i] : i.ToString();

    private static string Bar(double v, double max, int width)
    {
        double units = (max > 1e-300 ? Math.Clamp(v / max, 0.0, 1.0) : 0.0) * width;
        int full = Math.Min((int)Math.Floor(units), width);
        var sb = new StringBuilder(width);
        sb.Append('█', full);
        if (full < width)
        {
            int e = (int)Math.Round((units - full) * 8);
            if (e >= 8) sb.Append('█');
            else if (e > 0) sb.Append(BarEighths[e]);
        }
        return sb.ToString().PadRight(width);
    }

    private static char Glyph(double v, double max)
    {
        if (max <= 1e-300) return HeatRamp[0];
        int k = (int)Math.Round(Math.Clamp(v / max, 0.0, 1.0) * (HeatRamp.Length - 1));
        return HeatRamp[k];
    }

    private static char Cell(double fill)
    {
        if (fill >= 1.0) return '█';
        if (fill <= 0.0) return ' ';
        int e = (int)Math.Round(fill * 8);
        return e <= 0 ? ' ' : e >= 8 ? '█' : ColEighths[e - 1];
    }
}
