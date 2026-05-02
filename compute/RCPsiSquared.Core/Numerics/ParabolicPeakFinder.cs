namespace RCPsiSquared.Core.Numerics;

/// <summary>Peak-Q + sub-grid HWHM extraction for sampled curves.
///
/// <para>Locates the grid maximum, then refines the peak Q via a 3-point parabolic fit
/// y = a·x² + b·x + c (Q* = −b/(2a)). Half-maximum points on either side are linearly
/// interpolated between the bracketing grid points where K crosses K_max/2.</para>
///
/// <para>If the parabolic fit is concave-up (a ≥ 0) or extrapolates beyond the grid window,
/// the routine falls back to the grid maximum without sub-grid refinement.</para>
/// </summary>
public static class ParabolicPeakFinder
{
    public readonly record struct PeakInfo(double QPeak, double KMax, double? HwhmLeft, double? HwhmRight);

    public static PeakInfo Find(IReadOnlyList<double> qGrid, IReadOnlyList<double> kCurve)
    {
        if (qGrid.Count != kCurve.Count) throw new ArgumentException("qGrid and kCurve must have the same length.");
        if (qGrid.Count < 3) throw new ArgumentException("need at least 3 grid points.");

        int iMax = 0;
        double kMax = kCurve[0];
        for (int i = 1; i < kCurve.Count; i++)
        {
            if (kCurve[i] > kMax) { kMax = kCurve[i]; iMax = i; }
        }

        double qPeak = qGrid[iMax];

        // Parabolic refinement around the grid maximum, if interior.
        if (iMax > 0 && iMax < qGrid.Count - 1)
        {
            double x0 = qGrid[iMax - 1], x1 = qGrid[iMax], x2 = qGrid[iMax + 1];
            double y0 = kCurve[iMax - 1], y1 = kCurve[iMax], y2 = kCurve[iMax + 1];
            double[] coefs = FitQuadratic(x0, x1, x2, y0, y1, y2);
            double a = coefs[0], b = coefs[1], c = coefs[2];
            if (a < 0)
            {
                double qStar = -b / (2.0 * a);
                double kStar = c - b * b / (4.0 * a);
                if (Math.Abs(qStar - qGrid[iMax]) <= (qGrid[iMax + 1] - qGrid[iMax - 1]))
                {
                    qPeak = qStar;
                    kMax = kStar;
                }
            }
        }

        double half = kMax / 2.0;
        double? hwhmLeft = null;
        for (int i = iMax; i > 0; i--)
        {
            if (kCurve[i] < half)
            {
                double x0 = qGrid[i], x1 = qGrid[i + 1];
                double y0 = kCurve[i], y1 = kCurve[i + 1];
                double xHalf = x0 + (half - y0) * (x1 - x0) / (y1 - y0);
                hwhmLeft = qPeak - xHalf;
                break;
            }
        }
        double? hwhmRight = null;
        for (int i = iMax; i < qGrid.Count; i++)
        {
            if (kCurve[i] < half)
            {
                double x0 = qGrid[i - 1], x1 = qGrid[i];
                double y0 = kCurve[i - 1], y1 = kCurve[i];
                double xHalf = x0 + (half - y0) * (x1 - x0) / (y1 - y0);
                hwhmRight = xHalf - qPeak;
                break;
            }
        }
        return new PeakInfo(qPeak, kMax, hwhmLeft, hwhmRight);
    }

    /// <summary>Quadratic y = a·x² + b·x + c through three points. Returns [a, b, c].</summary>
    private static double[] FitQuadratic(double x0, double x1, double x2, double y0, double y1, double y2)
    {
        // Solve the 3x3 system by Cramer's rule.
        double det = (x0 - x1) * (x0 - x2) * (x1 - x2);
        double a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / det;
        double b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / det;
        double c = (x1 * x2 * (x1 - x2) * y0 + x2 * x0 * (x2 - x0) * y1 + x0 * x1 * (x0 - x1) * y2) / det;
        return new[] { a, b, c };
    }
}
