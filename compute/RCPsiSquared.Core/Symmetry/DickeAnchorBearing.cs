namespace RCPsiSquared.Core.Symmetry;

/// <summary>The role a <see cref="Knowledge.Claim"/> plays with respect to
/// the <see cref="DickeAnchor"/> 3-anchor set {Mirror, KIntermediate, Generic}.
/// Mirrors <see cref="F99AnchorRole"/> in structure (Direct / Parent /
/// Covers) but specialised to the symmetric Dicke superposition family at
/// γ ∈ {1, 1/2, 0} and α_total ∈ {0, 3/8, 1/2}.
///
/// <para>The DickeAnchor 3-set is a structural SUBSET of the F99 5-set
/// {0, 1/8, 1/4, 3/8, 1/2}: the three uniform-amplitude cases (c² = ∞, 1, 0
/// respectively) sit at three of the five F99 canonical anchors. The
/// remaining two F99 anchors (1/8 at γ=√3/2, 1/4 at γ=√2/2) correspond to
/// non-uniform Dicke c² ∈ {2√3+3, 1+√2} and are F99-only, not DickeAnchor.</para>
/// </summary>
public enum DickeAnchorRole
{
    /// <summary>Claim content IS a statement about one or more DickeAnchor
    /// enum cases directly (e.g. F98 asserts the K-intermediate long-time
    /// asymptote at uniform Dicke γ=1/2).</summary>
    Direct,

    /// <summary>Structural parent feeding the DickeAnchor inheritance graph
    /// without itself being a Dicke-specific statement.</summary>
    Parent,

    /// <summary>Claim covers a SET of DickeAnchor cases as a single theorem
    /// (e.g. F99 covers all three uniform-Dicke cases as a subset of its
    /// five canonical anchors).</summary>
    Covers,
}

/// <summary>Marker interface for <see cref="Knowledge.Claim"/>s that
/// participate in the <see cref="DickeAnchor"/> inheritance graph. Aggregated
/// by <see cref="DickeAnchorMap"/> into a 3-row coverage table {Mirror,
/// KIntermediate, Generic}, with gap-anchor detection analogous to
/// <see cref="F99AnchorMap"/>.
///
/// <para>The DickeAnchor map is a finer-grained reading of the same
/// territory: F99 covers 5 anchors via the full sin²(θ)/2 mechanism;
/// DickeAnchor covers the uniform-c subset of 3 anchors directly via the
/// X⊗N-eigenbasis decomposition. Both maps surface their gaps independently
/// — a claim can be a DickeAnchor-Direct without being F99-Direct (and vice
/// versa).</para>
/// </summary>
public interface IDickeAnchorBearing
{
    /// <summary>How this claim relates to the DickeAnchor 3-anchor set.</summary>
    DickeAnchorRole DickeRole { get; }

    /// <summary>The <see cref="DickeAnchor"/> cases this claim covers. Empty
    /// for <see cref="DickeAnchorRole.Parent"/>; singleton for the simplest
    /// <see cref="DickeAnchorRole.Direct"/>; multi-element for
    /// <see cref="DickeAnchorRole.Covers"/>.</summary>
    IReadOnlyList<DickeAnchor> DickeAnchors { get; }
}
