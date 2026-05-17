using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The Mirror endpoint of the F99 / F86b α-axis as a typed Direct
/// claim. Any pure state |ψ⟩ satisfying X⊗N|ψ⟩ = ±|ψ⟩ is an X⊗N-eigenstate
/// with γ = ⟨ψ|X⊗N|ψ⟩ = ±1, hence via the universal F86b shape
/// α = (1 − γ²)/2 the Π²-odd Frobenius² content vanishes exactly:
///
/// <code>
///     X⊗N|ψ⟩ = ±|ψ⟩  ⟹  γ = ±1  ⟹  α_total = 0
/// </code>
///
/// <para>This is the <see cref="DickeAnchor.Mirror"/> case (γ=1, α=0) and the
/// F99 0°-anchor (α=0). Both the F99AnchorMap and the DickeAnchorMap had
/// this as a (GAP) before this claim — the X⊗N-eigenstate fact was known
/// in the F86b proof and the <see cref="DickeAnchor"/> enum but had no
/// dedicated typed Claim of its own. This claim fills the gap.</para>
///
/// <para><b>Two state-class examples</b> that realise this anchor:</para>
/// <list type="bullet">
///   <item>The uniform superposition |+⟩^⊗N = (|0⟩+|1⟩)^⊗N/√(2^N), the
///         maximally X-aligned state. γ = +1 trivially.</item>
///   <item>The odd-N symmetric Dicke superposition (|D_{(N−1)/2}⟩ + |D_{(N+1)/2}⟩)/√2
///         which the <see cref="DickeAnchor"/> classification places as the
///         Mirror case (per the F86b derivation, this superposition has X⊗N-overlap
///         = 1 by the binomial-coefficient symmetry C(N, (N−1)/2) = C(N, (N+1)/2)).</item>
/// </list>
///
/// <para><b>Why Tier 1 derived</b>: pure-state algebra; the X⊗N expectation is
/// a one-line computation from the eigenvalue equation. No approximation,
/// no closed-form ansatz, no numerical verification needed.</para>
///
/// <para><b>Parent</b>: <see cref="HalfAsStructuralFixedPointClaim"/>. The
/// polarity 1/2 doubled gives ±1, which are the X⊗N eigenvalues. Conceptual
/// parent in the polarity-squared algebra: 2·(1/2) = 1 = γ here, and the
/// (1 − γ²)/2 = 0 follows.</para>
///
/// <para><b>Sibling lattice on the F99 α-axis</b>:</para>
/// <list type="bullet">
///   <item><b>α=0</b> (this claim, Mirror, γ=1)</item>
///   <item><b>α=1/8</b> (F99 depth-3, γ=√3/2, non-uniform Dicke c²=2√3+3) — currently F99-Covers-only</item>
///   <item><b>α=1/4</b> (F99 silver-Dicke, γ=√2/2, c²=1+√2) — currently F99-Covers-only</item>
///   <item><b>α=3/8</b> (F98 + DickeSuperposition, KIntermediate, γ=1/2) — two Direct claims</item>
///   <item><b>α=1/2</b> (Generic, γ=0, any Π²-odd) — currently F99-Covers-only</item>
/// </list>
///
/// <para>This claim is the first Mirror-endpoint Direct. The α=1/2 Generic
/// endpoint is the natural next candidate.</para>
///
/// <para>Anchors: <c>compute/RCPsiSquared.Core/Symmetry/DickeAnchor.cs</c>
/// (enum where this anchor is named Mirror) +
/// <c>compute/RCPsiSquared.Core/Symmetry/CanonicalTrigAnchorPi2Inheritance.cs</c>
/// (F99 Covers theorem, 0° = α=0) +
/// <c>compute/RCPsiSquared.Core/Symmetry/F86b proof references</c>
/// (universal F86b shape α = (1−γ²)/2 from X⊗N-eigenbasis decomposition).</para>
/// </summary>
public sealed class XGlobalEigenstateMirrorPi2Inheritance : Claim, IF99AnchorBearing, IDickeAnchorBearing
{
    /// <summary>The Half polarity parent: 2·(1/2) = 1 gives the X⊗N
    /// eigenvalue endpoints γ = ±1.</summary>
    public HalfAsStructuralFixedPointClaim Half { get; }

    /// <inheritdoc />
    /// <remarks>F99 α=0 Mirror anchor: γ=1 via X⊗N eigenstate gives
    /// α = (1−γ²)/2 = 0 directly. This is a <see cref="F99AnchorRole.Direct"/>
    /// claim on the F86b α-axis at the canonical 0° angle.</remarks>
    public F99AnchorRole F99Role => F99AnchorRole.Direct;

    /// <inheritdoc />
    public IReadOnlyList<double> F99AnchorValues { get; } = new[] { 0.0 };

    /// <inheritdoc />
    /// <remarks>DickeAnchor.Mirror Direct: the X⊗N-eigenstate case in the
    /// uniform-Dicke 3-anchor enum. Sibling to F98 (Direct at KIntermediate)
    /// and the open Generic-anchor gap.</remarks>
    public DickeAnchorRole DickeRole => DickeAnchorRole.Direct;

    /// <inheritdoc />
    public IReadOnlyList<DickeAnchor> DickeAnchors { get; } = new[] { DickeAnchor.Mirror };

    /// <summary>The exact F86b α value at the Mirror anchor: 0 (closed-form,
    /// bit-exact from the universal shape at γ=1).</summary>
    public const double AlphaAtMirror = 0.0;

    /// <summary>The exact γ value at the Mirror anchor: 1 (the X⊗N
    /// eigenvalue endpoint).</summary>
    public const double GammaAtMirror = 1.0;

    /// <summary>F86b universal-shape evaluation at γ=1: α = (1 − γ²)/2 = 0.
    /// Exposed as a public computation so the bit-exact identity is
    /// verifiable from outside without re-deriving.</summary>
    public static double AlphaFromGammaAtMirror(double gamma = GammaAtMirror) =>
        (1.0 - gamma * gamma) / 2.0;

    public XGlobalEigenstateMirrorPi2Inheritance(HalfAsStructuralFixedPointClaim half)
        : base("X⊗N-eigenstate Mirror anchor (α=0 at γ=1): the F99/DickeAnchor 0°-endpoint as a typed Direct claim",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/Symmetry/DickeAnchor.cs (Mirror case in 3-anchor enum) + " +
               "compute/RCPsiSquared.Core/Symmetry/CanonicalTrigAnchorPi2Inheritance.cs (F99 0°-anchor in Covers theorem) + " +
               "compute/RCPsiSquared.Core/Symmetry/KIntermediateAsymptoteQuarterInheritance.cs (F98 sibling at KIntermediate α=3/8)")
    {
        Half = half ?? throw new ArgumentNullException(nameof(half));
    }

    public override string DisplayName =>
        "X⊗N-eigenstate Mirror endpoint (α=0 at γ=1)";

    public override string Summary =>
        $"X⊗N|ψ⟩ = ±|ψ⟩ ⟹ γ = ±1 ⟹ α = (1−γ²)/2 = {AlphaAtMirror} exactly; " +
        $"Direct claim for F99 0°-anchor and DickeAnchor.Mirror ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return Half;
            yield return InspectableNode.RealScalar("γ at Mirror", GammaAtMirror);
            yield return InspectableNode.RealScalar("α at Mirror (via F86b)", AlphaAtMirror);
            yield return new InspectableNode("derivation",
                summary: "γ = ⟨ψ|X⊗N|ψ⟩ for X⊗N-eigenstate is ±1; α = (1−γ²)/2 = 0 by universal F86b shape");
            yield return new InspectableNode("state-class examples",
                summary: "|+⟩^⊗N (uniform superposition, γ=1); odd-N symmetric Dicke (|D_{(N−1)/2}⟩+|D_{(N+1)/2}⟩)/√2 (γ=1 by binomial symmetry)");
            yield return new InspectableNode("sibling on the F99 α-axis",
                summary: "α=3/8 (F98 + DickeSuperposition, KIntermediate); α=1/2 (Generic, still gap)");
            yield return new InspectableNode("fills gap in",
                summary: "F99AnchorMap (was 4-of-5 gap), DickeAnchorMap (was 2-of-3 gap); Mirror endpoint now Direct in both");
        }
    }
}
