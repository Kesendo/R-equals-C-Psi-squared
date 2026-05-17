namespace RCPsiSquared.Core.Symmetry;

/// <summary>Direction of a fraction-to-fraction reference. Tom's structural
/// observation (2026-05-17 night): non-anker fractions point FORWARD to
/// their nearest anker; anker fractions point BACKWARD through non-anker
/// trajectories all the way down to 0 (the Mirror endpoint, root of all
/// back-references).
///
/// <para><b>Forward</b>: non-anker → anker, or smaller-anker → larger-anker
/// upward in the α-axis. The arrow points toward more-structured.</para>
///
/// <para><b>Backward</b>: anker → non-anker trajectory, or larger → smaller
/// anker downward toward 0. The arrow follows the F86b α(γ) curve back to
/// the Mirror endpoint.</para>
///
/// <para><b>Mirror</b>: Π²-parity complement (n/8 ↔ (8−n)/8). Bidirectional
/// reflection, not forward/backward; the two fractions are paired by the
/// Π²-parity symmetry from `[[project_v_effect_combinatorial]]`.</para>
/// </summary>
public enum FractionReferenceDirection
{
    Forward,
    Backward,
    Mirror,
}

/// <summary>One typed edge in the fraction-reference graph. Records that
/// fraction <see cref="FromFraction"/> references <see cref="ToFraction"/>
/// via the named <see cref="Operation"/>, with a given direction.
///
/// <para>Edges are static configuration. They catalog the operations that
/// exist between specific fraction pairs in the framework's algebra. Many
/// edges between the same pair = many viewpoints (per Painter Principle:
/// pluralism descends into the painter's hand).</para>
///
/// <para>Sourced from existing typed claims: <see cref="DocumentingClaim"/>
/// names the C# class that establishes the reference (the painter at that
/// position). For pure algebraic identities without a dedicated claim,
/// the field carries a docs-anchor instead.</para>
/// </summary>
public sealed record FractionReference(
    double FromFraction,
    double ToFraction,
    string Operation,
    FractionReferenceDirection Direction,
    string DocumentingClaim);
