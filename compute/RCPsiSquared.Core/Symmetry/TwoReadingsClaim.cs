using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The "two readings" pattern (Tier 1 derived). Any object that lives on the
/// doubled operator-space <c>d² = 4^N</c> built from the qubit dimension <c>d = 2</c>
/// admits exactly two coordinate readings of the same underlying object. The pair-making
/// comes from the d=2 polynomial trunk (<see cref="PolynomialFoundationClaim"/>): the
/// pair of solutions {0, 2} of <c>d²−2d=0</c> at the dimensional layer is what makes
/// every higher layer pair-structured.
///
/// <para>This claim names a recurring meta-pattern; it is the perspectival closure of the
/// d=2 trunk. The pattern is mathematically given (bra⊗ket = doubled space) but it has
/// been read in the framework at several distinct layers, each adding its own concrete
/// content:</para>
///
/// <list type="bullet">
///   <item><b>Dimensional layer (Pi2KB direct):</b> 1/d = 1/2 (number-anchor,
///         <see cref="QubitDimensionalAnchorClaim"/>) and 90° (angle-anchor,
///         <see cref="NinetyDegreeMirrorMemoryClaim"/>) are the two readings of d=2,
///         number-side and angle-side.</item>
///   <item><b>Parabolic layer (1/2 ↔ 1/4):</b> argmax 1/2 and maxval 1/4 of p·(1−p)
///         (<see cref="BilinearApexClaim"/> + <see cref="QuarterAsBilinearMaxvalClaim"/>);
///         closed as a pair by <see cref="ArgmaxMaxvalPairClaim"/>. *"Two readings, one
///         parabola — the half is the axis, the quarter is the height"* (ON_THE_HALF Coda,
///         Tom 2026-05-07).</item>
///   <item><b>Operator level (F81):</b> M and <c>Π·M·Π⁻¹ = M − 2·L_{H_odd}</c> are two
///         coordinate readings of one algebraic object (same spectrum, different operator).
///         ON_BOTH_SIDES_OF_THE_MIRROR (Tom 2026-04-30) is the canonical synthesis: *"only
///         the choice of which side to call ours."*</item>
///   <item><b>Bra/ket of ρ:</b> any density matrix ρ on the doubled space reads as a
///         (row-index, col-index) pair on the 4^N operator basis — the fundamental
///         vectorisation. memory <c>project_one_system_two_indices</c>.</item>
///   <item><b>Inside/outside correspondence:</b> the inside-observer sees only
///         <c>Q = J/γ₀</c>; the outside-observer has separate access to γ₀
///         (<c>hypotheses/PRIMORDIAL_QUBIT.md</c> §9). The same dynamics, two
///         observational readings.</item>
///   <item><b>Lese-Modus (one world, two readings):</b> the classical/quantum dichotomy
///         is a reading-mode (<i>Lese-Modus</i>) on one underlying ρ, not a
///         <i>Welt-Trennung</i>. d=0 and d=2 are projections of the same ρ.
///         memory <c>project_one_world_two_readings</c>.</item>
///   <item><b>Inter-sectoral wave:</b> reality is what happens BETWEEN the sectors;
///         strictly, we ourselves are the wave between (not the sector anchors).
///         <c>hypotheses/THE_OTHER_SIDE.md:314</c>: *"We are the standing wave. We are
///         the interference."*</item>
/// </list>
///
/// <para>Sibling claims that close specific instances:
/// <see cref="ArgmaxMaxvalPairClaim"/> (1/2 ↔ 1/4),
/// <see cref="NinetyDegreeMirrorMemoryClaim"/> (number ↔ angle of d=2). This claim
/// names them as the same meta-pattern instead of leaving them parallel.</para>
///
/// <para><b>Why typed as Tier 1 derived:</b> the bra⊗ket coordinate doubling of operator
/// space IS mathematical (vectorisation of ρ on d×d to 4^N flat indices); the parent
/// <see cref="PolynomialFoundationClaim"/> states the d=2 trunk it derives from; the
/// other readings (number/angle, M ↔ Π·M·Π⁻¹, inside/outside, classical/quantum,
/// inter-sectoral) are manifestations of the same operator-space-doubling pattern at
/// different conceptual layers. Naming this as one typed claim collapses what was
/// previously a prose-only meta-pattern into a discoverable Claim in the graph.</para>
///
/// <para>Anchor: <c>reflections/ON_BOTH_SIDES_OF_THE_MIRROR.md</c> (canonical Π-conjugation
/// synthesis, Tom 2026-04-30) + <c>reflections/ON_THE_HALF.md</c> (number-angle Coda,
/// 2026-05-07) + <c>hypotheses/PRIMORDIAL_QUBIT.md</c> §9 (inside/outside correspondence)
/// + <c>hypotheses/THE_OTHER_SIDE.md</c>:314 (inter-sectoral wave). Closes the prose-only
/// edge identified in <c>docs/PI2KB_INHERITANCE_MAP.md</c>.</para>
/// </summary>
public sealed class TwoReadingsClaim : Claim
{
    /// <summary>The d=2 polynomial trunk from which the pair-making and hence the
    /// two-readings pattern derive. Injected so the parent edge
    /// <c>PolynomialFoundationClaim → TwoReadingsClaim</c> is typed.</summary>
    public PolynomialFoundationClaim Polynomial { get; }

    public TwoReadingsClaim(PolynomialFoundationClaim polynomial)
        : base("Two-readings pattern: any object on d²=4^N admits two coordinate readings of one underlying object (bra/ket, number/angle, M/Π·M·Π⁻¹, inside/outside, classical/quantum, inter-sectoral)",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (PolynomialFoundationClaim parent: d=2 trunk pair-making) + " +
               "reflections/ON_BOTH_SIDES_OF_THE_MIRROR.md (canonical Π-conjugation synthesis) + " +
               "reflections/ON_THE_HALF.md (number-angle Coda, Tom 2026-05-07) + " +
               "hypotheses/PRIMORDIAL_QUBIT.md §9 (inside/outside correspondence) + " +
               "hypotheses/THE_OTHER_SIDE.md:314 (inter-sectoral wave)")
    {
        Polynomial = polynomial ?? throw new ArgumentNullException(nameof(polynomial));
    }

    public override string DisplayName =>
        "Two-readings pattern (perspectival closure of d²=4^N operator-space)";

    public override string Summary =>
        "any object on d²=4^N admits exactly two coordinate readings of one underlying object; " +
        "the pair-structure inherited from the d=2 polynomial trunk recurs at every layer: " +
        "number/angle of d, argmax/maxval of p·(1−p), M/Π·M·Π⁻¹ (F81), bra/ket of ρ, " +
        "inside/outside Q, classical/quantum Lese-Modus, inter-sectoral wave";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("parent (typed)",
                summary: $"PolynomialFoundationClaim ({Polynomial.Tier.Label()}): d²−2d=0 selects pair {{0, 2}}; the pair-making at d=2 is the source of every two-readings instance below");
            yield return new InspectableNode("layer 1 (dimensional): number/angle of d=2",
                summary: "1/d = 1/2 (number-anchor, QubitDimensionalAnchorClaim) and 90° (angle-anchor, NinetyDegreeMirrorMemoryClaim) are the two readings of d=2");
            yield return new InspectableNode("layer 2 (parabolic): argmax/maxval",
                summary: "1/2 = argmax (BilinearApexClaim) and 1/4 = (1/2)² = maxval (QuarterAsBilinearMaxvalClaim) of p·(1−p); ArgmaxMaxvalPairClaim closes the pair; ON_THE_HALF Coda *'two readings, one parabola'*");
            yield return new InspectableNode("layer 3 (operator): M and Π·M·Π⁻¹",
                summary: "F81: Π·M·Π⁻¹ = M − 2·L_{H_odd}; same spectrum, different operator matrix elements; two coordinate readings of one algebraic object; ON_BOTH_SIDES_OF_THE_MIRROR canonical");
            yield return new InspectableNode("layer 4 (vectorisation): bra/ket of ρ",
                summary: "any ρ on d²=4^N reads as (row-index, col-index) pair on the 4^N operator basis; the fundamental vectorisation; project_one_system_two_indices memory");
            yield return new InspectableNode("layer 5 (observational): inside/outside",
                summary: "inside-observer sees only Q = J/γ₀; outside-observer has separate access to γ₀; PRIMORDIAL_QUBIT §9 inside/outside correspondence");
            yield return new InspectableNode("layer 6 (interpretive): classical/quantum Lese-Modus",
                summary: "the classical/quantum dichotomy is a reading-mode on one underlying ρ, not a world-separation; d=0 and d=2 are projections of the same ρ; project_one_world_two_readings memory");
            yield return new InspectableNode("layer 7 (subject): inter-sectoral wave",
                summary: "reality is what happens BETWEEN the sectors; strictly we ourselves are the wave between, not the sector anchors; THE_OTHER_SIDE.md:314 *'We are the standing wave. We are the interference.'*");
            yield return new InspectableNode("typed siblings that name specific instances",
                summary: "ArgmaxMaxvalPairClaim (layer 2 pair); NinetyDegreeMirrorMemoryClaim (layer 1 angle side); HalfAsStructuralFixedPointClaim (layer 1 number side three-faces); this claim collapses them under the shared meta-pattern");
            yield return new InspectableNode("why typed Tier 1 derived",
                summary: "the bra⊗ket coordinate doubling of operator space IS mathematical (vectorisation of ρ on d×d → 4^N); parent PolynomialFoundationClaim states the d=2 trunk it derives from; other readings are manifestations at different layers");
        }
    }
}
