using System.Globalization;

namespace RCPsiSquared.Core.Numerics;

/// <summary>Formatting helpers used by the F86 inspection layer to produce signed-delta
/// strings (Unicode minus + ASCII plus) without requiring a Custom-Format-String section
/// dance that produces "−+0.0000" for negative zeros.
/// </summary>
public static class Formatting
{
    /// <summary>Format a value as a signed delta with a Unicode minus for negatives and an
    /// ASCII plus for non-negatives, e.g. "+0.0050" or "−0.0105".</summary>
    public static string SignedDelta(double value, string format = "F4") =>
        (value >= 0 ? "+" : "−") + Math.Abs(value).ToString(format, CultureInfo.CurrentCulture);
}
