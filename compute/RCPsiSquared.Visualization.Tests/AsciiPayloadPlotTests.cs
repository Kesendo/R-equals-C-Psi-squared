using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Visualization.Inspection;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Visualization.Tests;

/// <summary>Behaviour of the ASCII payload renderers that turn an <see cref="InspectablePayload"/>
/// into terminal art for the <see cref="ConsoleTreeRenderer"/>: Vector -> magnitude bars,
/// MatrixView -> magnitude heatmap, Curve -> multi-row bar plot. The assertions check structure
/// (bar length ordering, full/blank glyphs, downsampling, axis labels), not exact strings, so the
/// glyph ramps can be tuned without churning the tests.</summary>
public class AsciiPayloadPlotTests
{
    private static ComplexVector Vec(params double[] re) =>
        ComplexVector.Build.DenseOfArray(re.Select(r => new Complex(r, 0)).ToArray());

    // ---- Vector: magnitude bars --------------------------------------------------------------

    [Fact]
    public void Vector_BarLengthScalesWithMagnitude()
    {
        var lines = AsciiPayloadPlot.Vector(Vec(1.0, 0.0, 0.5));
        Assert.Equal(3, lines.Count);
        int full0 = lines[0].Count(ch => ch == '█');
        int full1 = lines[1].Count(ch => ch == '█');
        int full2 = lines[2].Count(ch => ch == '█');
        Assert.True(full0 > full2, "max-magnitude row should have the longest bar");
        Assert.True(full2 > full1, "half-magnitude row longer than the zero row");
        Assert.Equal(0, full1);
    }

    [Fact]
    public void Vector_UsesComponentLabels()
    {
        var lines = AsciiPayloadPlot.Vector(Vec(1.0, 1.0), new[] { "qA", "qB" });
        Assert.StartsWith("qA", lines[0]);
        Assert.StartsWith("qB", lines[1]);
    }

    [Fact]
    public void Vector_AllZero_NoBars()
    {
        var lines = AsciiPayloadPlot.Vector(Vec(0.0, 0.0, 0.0));
        Assert.All(lines, l => Assert.DoesNotContain('█', l));
    }

    [Fact]
    public void Vector_LongVector_TruncatesToTopByMagnitude()
    {
        var big = Vec(Enumerable.Range(0, 40).Select(i => (double)i).ToArray());
        var lines = AsciiPayloadPlot.Vector(big, null, width: 64, maxRows: 8);
        Assert.Equal(9, lines.Count); // 8 rows + the "(top 8 of 40 ...)" note
        Assert.Contains("of 40", lines[^1]);
    }

    // ---- Matrix: magnitude heatmap -----------------------------------------------------------

    [Fact]
    public void Matrix_MaxCellRendersFullBlock_WithMaxFooter()
    {
        var m = ComplexMatrix.Build.DenseOfArray(new Complex[,] { { 2, 0 }, { 0, 1 } });
        var lines = AsciiPayloadPlot.Matrix(m);
        Assert.Contains(lines, l => l.Contains('█'));        // the max-magnitude cell (2)
        Assert.Contains("max 2", string.Join("\n", lines));  // footer carries the scale
    }

    [Fact]
    public void Matrix_ZeroMatrix_NoFullBlock()
    {
        var m = ComplexMatrix.Build.Dense(3, 3, Complex.Zero);
        var lines = AsciiPayloadPlot.Matrix(m);
        Assert.DoesNotContain('█', string.Join("", lines));
    }

    [Fact]
    public void Matrix_LargeMatrix_DownsampledToCap()
    {
        var m = ComplexMatrix.Build.Dense(40, 40, Complex.One);
        var lines = AsciiPayloadPlot.Matrix(m, width: 64, maxDim: 24);
        Assert.Equal(25, lines.Count); // 24 grid rows + footer
        Assert.Contains("40×40", string.Join("\n", lines));
    }

    // ---- Curve: multi-row bar plot -----------------------------------------------------------

    [Fact]
    public void Curve_BottomRowFillsMoreThanTop_ForMonotonic()
    {
        var y = new double[] { 0, 1, 2, 3, 4 };
        var x = new double[] { 0, 1, 2, 3, 4 };
        var lines = AsciiPayloadPlot.Curve(x, y, "Q", width: 48, height: 4);
        Assert.Equal(5, lines.Count); // 4 rows + x-axis line
        int bottomFull = lines[3].Count(ch => ch == '█'); // row = 1, the bottom data row
        int topFull = lines[0].Count(ch => ch == '█');     // row = height, the top
        Assert.True(bottomFull > topFull, "a rising curve fills more at the bottom than the top");
    }

    [Fact]
    public void Curve_HasYAxisLabelsAndXRange()
    {
        var y = new double[] { 0.01, 0.2, 0.42, 0.2, 0.01 };
        var x = new double[] { 0.2, 1.0, 1.5, 3.0, 4.0 };
        var lines = AsciiPayloadPlot.Curve(x, y, "Q", height: 4);
        Assert.Contains("0.42", lines[0]);   // ymax on the top row
        Assert.Contains("0.01", lines[3]);   // ymin on the bottom data row
        Assert.Contains("Q", lines[^1]);     // x-axis label
        Assert.Contains("4", lines[^1]);     // x range end
    }

    // ---- dispatch ----------------------------------------------------------------------------

    [Fact]
    public void Render_EmptyForScalarPayloads_NonEmptyForVector()
    {
        Assert.Empty(AsciiPayloadPlot.Render(new InspectablePayload.Real("r", 1.0)));
        Assert.Empty(AsciiPayloadPlot.Render(InspectablePayload.Empty));
        Assert.NotEmpty(AsciiPayloadPlot.Render(new InspectablePayload.Vector("v", Vec(1.0, 0.5))));
    }
}
