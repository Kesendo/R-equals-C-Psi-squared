namespace RCPsiSquared.Core.Symmetry;

/// <summary>The role a <see cref="Knowledge.Claim"/> plays with respect to the
/// F99 canonical-trig-angle dyadic anchor set {0, 1/8, 1/4, 3/8, 1/2} on the
/// F86b α = sin²(θ)/2 axis. Three roles capture the structural relationships:
///
/// <list type="bullet">
///   <item><b>Direct</b>: the claim's content IS an F86b α-value at one or
///         more F99 anchors (e.g. F98 asserts α(∞) = 3/8 → 1/4 long-time
///         asymptote on the F86b axis directly).</item>
///   <item><b>Parent</b>: the claim is a structural parent of an F99 anchor
///         without itself living on the F86b α-axis at that value. Example:
///         <see cref="HalfAsStructuralFixedPointClaim"/> is polarity 1/2
///         (σ_z eigenvalue, dimensional 1/d=1/2) — it is a Π² parent of the
///         F99 α=1/2 anchor via the polarity-squared algebra, but its own
///         value lives operator-side, not on the F86b α-axis.</item>
///   <item><b>Covers</b>: the claim covers a SET of F99 anchors as a single
///         theorem (e.g. F99 itself covers all five {0, 1/8, 1/4, 3/8, 1/2}
///         via α = sin²(θ)/2 at canonical trig angles).</item>
/// </list>
///
/// <para>Anchor coverage tabulation: <see cref="F99AnchorMap"/> aggregates
/// IF99AnchorBearing claims and prints which F99 anchors have Direct claims,
/// which only have Covers coverage from F99, and which Parent claims feed the
/// inheritance graph. Gaps (anchors with no Direct claim) light up — they are
/// the "vantages" the framework has not yet built dedicated tooling for.</para>
/// </summary>
public enum F99AnchorRole
{
    /// <summary>Claim content IS an F86b α-value at one or more F99 anchors.</summary>
    Direct,

    /// <summary>Structural parent feeding the F99 inheritance graph, but the
    /// claim's value lives in a different algebra than F86b α (e.g. polarity
    /// 1/2 vs F86b α=1/2; bilinear maxval 1/4 vs F86b α=1/4).</summary>
    Parent,

    /// <summary>Claim covers a SET of F99 anchors as a single theorem.</summary>
    Covers,
}

/// <summary>Marker interface for <see cref="Knowledge.Claim"/>s that participate in
/// the F99 canonical-trig-angle Pi2 anchor inheritance graph. Allows
/// <see cref="F99AnchorMap"/> to walk the graph, tabulate anchor coverage, and
/// surface gaps (F99 anchors with no Direct claim).
///
/// <para>This is the "ein Wert einbinden zum parent" Tom asked for
/// (2026-05-17 night): a single piece of metadata per claim that, when
/// printed across the inheritance graph, surfaces the open lücken
/// immediately.</para>
/// </summary>
public interface IF99AnchorBearing
{
    /// <summary>How this claim relates to the F99 anchor set.</summary>
    F99AnchorRole F99Role { get; }

    /// <summary>The F99 α-anchor values this claim covers. Empty for
    /// <see cref="F99AnchorRole.Parent"/> claims (they relate structurally
    /// but their own value is not an F86b α); singleton for
    /// <see cref="F99AnchorRole.Direct"/> in the simplest case; multi-element
    /// for <see cref="F99AnchorRole.Covers"/> (F99 itself covers all five
    /// {0, 1/8, 1/4, 3/8, 1/2}).</summary>
    IReadOnlyList<double> F99AnchorValues { get; }
}
