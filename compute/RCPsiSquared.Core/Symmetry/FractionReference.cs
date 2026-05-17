namespace RCPsiSquared.Core.Symmetry;

/// <summary>Direction of a fraction-to-fraction reference. Tom's structural
/// observation (2026-05-17 night): non-anker fractions point FORWARD to
/// their nearest anker; anker fractions point BACKWARD through non-anker
/// trajectories on the α-axis. The α-axis is the folded picture under
/// α = (1−γ²)/2; the true polarity-mirror structure lives on the γ-axis
/// (in a separate <c>PolarityMirrorMap</c>), and 0 is NOT a root but the
/// convergence point where the ±γ-polarity sides meet under folding (per
/// <see cref="PolarityLayerOriginClaim"/>).
///
/// <para><b>Forward</b>: non-anker → anker, or smaller-anker → larger-anker
/// upward in the α-axis. The arrow points toward more-structured.</para>
///
/// <para><b>Backward</b>: anker → non-anker trajectory, or larger → smaller
/// anker downward toward 0. The arrow follows the F86b α(γ) curve toward
/// the Mirror convergence point (where ±γ both give α=0).</para>
///
/// <para><b>Mirror</b>: Π²-parity complement (n/8 ↔ (8−n)/8). Bidirectional
/// reflection, not forward/backward; the two fractions are paired by the
/// Π²-parity symmetry from `[[project_v_effect_combinatorial]]`.</para>
///
/// <para><b>Polarity</b>: marks an edge where the α=value on this axis is
/// reached from BOTH ±γ-polarity sides under the (1−γ²)/2 folding. Self-loop
/// at α=0 (Mirror convergence: γ=±1) is the canonical case. Distinct from
/// Mirror (Π²-parity) which lives across n/8 complements on the α-axis.</para>
/// </summary>
public enum FractionReferenceDirection
{
    Forward,
    Backward,
    Mirror,
    Polarity,
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
