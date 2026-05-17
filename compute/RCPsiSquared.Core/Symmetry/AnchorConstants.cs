using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Shared constants and helpers for the F99 / Dicke / polarity
/// anchor map family. Centralises the small primitives that were scattered
/// across <see cref="F99AnchorMap"/>, <see cref="DickeAnchorMap"/>,
/// <see cref="FractionReferenceGraph"/>, and <see cref="PolarityMirrorMap"/>:
///
/// <list type="bullet">
///   <item><see cref="Tol"/>: the canonical floating-point tolerance for
///         anchor-value equality checks across all maps.</item>
///   <item><see cref="ShortClaimName"/>: strips conventional suffixes
///         (<c>Pi2Inheritance</c>, <c>Claim</c>) from a Claim type name for
///         compact map rendering.</item>
///   <item><see cref="FormatEighthFraction"/>: pretty-prints α-axis values
///         on the dyadic ladder {0, 1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8, 1}.</item>
/// </list>
/// </summary>
public static class AnchorConstants
{
    /// <summary>Canonical floating-point tolerance for anchor equality.
    /// Used by every map that compares double-valued anchor positions.</summary>
    public const double Tol = 1e-12;

    /// <summary>Strip conventional Claim-naming suffixes for compact display
    /// in map render output. <c>F98KIntermediateAsymptoteQuarterInheritance</c>
    /// becomes <c>F98KIntermediateAsymptoteQuarter</c> (drops <c>Inheritance</c>);
    /// <c>HalfAsStructuralFixedPointClaim</c> becomes
    /// <c>HalfAsStructuralFixedPoint</c> (drops <c>Claim</c>).</summary>
    public static string ShortClaimName(Claim c)
    {
        var t = c.GetType().Name;
        if (t.EndsWith("Pi2Inheritance"))
            t = t[..^"Pi2Inheritance".Length];
        if (t.EndsWith("Claim"))
            t = t[..^"Claim".Length];
        return t;
    }

    /// <summary>Pretty-print common eighth-fractions; falls back to G6 numeric
    /// format for anything else. Covers the full dyadic-ladder positions
    /// surfaced by the F99 anker set plus its Π²-parity complements.</summary>
    public static string FormatEighthFraction(double value)
    {
        if (Math.Abs(value)             < Tol) return "0";
        if (Math.Abs(value - 1.0 / 8.0) < Tol) return "1/8";
        if (Math.Abs(value - 1.0 / 4.0) < Tol) return "1/4";
        if (Math.Abs(value - 3.0 / 8.0) < Tol) return "3/8";
        if (Math.Abs(value - 1.0 / 2.0) < Tol) return "1/2";
        if (Math.Abs(value - 5.0 / 8.0) < Tol) return "5/8";
        if (Math.Abs(value - 3.0 / 4.0) < Tol) return "3/4";
        if (Math.Abs(value - 7.0 / 8.0) < Tol) return "7/8";
        if (Math.Abs(value - 1.0)       < Tol) return "1";
        return value.ToString("G6");
    }
}
