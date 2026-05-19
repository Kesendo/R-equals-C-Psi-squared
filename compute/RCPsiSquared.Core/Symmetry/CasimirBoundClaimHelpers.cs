namespace RCPsiSquared.Core.Symmetry;

/// <summary>Small primitives shared by the bipartite-complete Casimir Im-max bound
/// sister claims (<see cref="RingN4DihedralLockClaim"/>, <see cref="StarImMaxBoundClaim"/>)
/// from the 2026-05-19 Q-sweep sprint:
///
/// <list type="bullet">
///   <item><see cref="QSweepAnchorLabels"/>: the canonical six (label, Q) pairs used in
///         the 2026-05-19 Q-sweep at γ₀=0.05 (Q ∈ {0.5, 1.0, 1.5, √3, 2.0, 2.5}). The same
///         six Q-values anchor both Ring N=4 and Star Im-max bound claims; the label
///         text is the project-wide convention (matches the Role strings in
///         <see cref="QAnchorMap.CanonicalAnchors"/> for the six anchors that appear in
///         the sweep, dropping the four non-sweep anchors 0.2, 0.35, 1.2, 1.6, 1.8).</item>
///   <item><see cref="QSqrt3"/>: <c>Math.Sqrt(3.0) ≈ 1.7321</c> as a typed double rather
///         than a hard-coded 16-digit literal. Matches
///         <see cref="LindbladAbsorptionMatchAtSixtyDegreesClaim.QValue"/>; named
///         locally so the Casimir bound claims do not pull a cross-family dependency
///         on F95 / Absorption Theorem just for the constant.</item>
///   <item><see cref="RequireFiniteNonNegative"/>: validation guard reused by the
///         <c>Predict</c> / <c>PredictImOverSigma</c> entry points on both claims.
///         Throws <see cref="ArgumentException"/> with a consistent message.</item>
/// </list>
/// </summary>
public static class CasimirBoundClaimHelpers
{
    /// <summary><c>Math.Sqrt(3.0)</c> as a shared constant so the two Casimir-bound
    /// claims do not each hard-code the 16-digit literal. Matches the
    /// <see cref="LindbladAbsorptionMatchAtSixtyDegreesClaim.QValue"/> definition.</summary>
    public static readonly double QSqrt3 = Math.Sqrt(3.0);

    /// <summary>The six canonical (Q label, Q value) pairs used by the 2026-05-19
    /// Q-sweep at γ₀=0.05. Same list anchors both
    /// <see cref="RingN4DihedralLockClaim.EmpiricalAnchors"/> (6 rows) and
    /// <see cref="StarImMaxBoundClaim.EmpiricalAnchors"/> (24 rows = 4 N × 6 Q).
    /// Labels follow the QAnchorMap convention.</summary>
    public static IReadOnlyList<(string Label, double Q)> QSweepAnchorLabels { get; } =
        new (string, double)[]
        {
            ("Q=0.5 sub-balance",    0.5),
            ("Q=1.0 Balance",        1.0),
            ("Q=1.5 F86 Q_peak c=2", 1.5),
            ("Q=√3 canonical 60°",   QSqrt3),
            ("Q=2.0 Q_EP idealized", 2.0),
            ("Q=2.5 Endpoint orbit", 2.5),
        };

    /// <summary>Throws <see cref="ArgumentException"/> if <paramref name="value"/> is
    /// not finite or is negative. Used by the Predict / PredictImOverSigma entry
    /// points on both Casimir-bound claims to validate <c>J</c> and <c>Q</c> inputs.</summary>
    public static void RequireFiniteNonNegative(double value, string paramName)
    {
        if (!double.IsFinite(value))
            throw new ArgumentException(
                $"{paramName} must be finite; got {value}.", paramName);
        if (value < 0)
            throw new ArgumentException(
                $"{paramName} must be non-negative; got {value}.", paramName);
    }
}
