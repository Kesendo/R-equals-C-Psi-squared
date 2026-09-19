using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>
/// A discoverable open synthesis for several useful pairs of descriptions in the
/// project. The exact examples keep their own proofs: vectorisation supplies row
/// and column coordinates, a parabola has distinct argmax and maxval features, and
/// F81 supplies a specific Π-conjugation identity. The observational and reflective
/// examples are explicitly interpretive comparisons.
///
/// <para>Living on a <c>d²</c>-dimensional operator space does not imply that an
/// arbitrary object has exactly two readings. Nor does the polynomial
/// <c>d²-2d=0</c> make the examples below descendants of one mathematical theorem.
/// The stable type remains in the registry so the possible synthesis can be asked
/// about without recording a false ancestry edge.</para>
/// </summary>
public sealed class TwoReadingsClaim : Claim
{
    public TwoReadingsClaim()
        : base(
            "Several exact constructions and interpretive comparisons come in pairs, " +
            "but d² operator-space vectorisation does not imply that every object has exactly two readings",
            Tier.OpenQuestion,
            "reflections/ON_BOTH_SIDES_OF_THE_MIRROR.md (specific Π-conjugation synthesis) + " +
            "reflections/ON_THE_HALF.md (argmax/maxval comparison) + " +
            "hypotheses/PRIMORDIAL_QUBIT.md §9 and hypotheses/THE_OTHER_SIDE.md (interpretive catalogue)")
    {
    }

    public override string DisplayName =>
        "Two readings? (open synthesis, not a universal d² consequence)";

    public override string Summary =>
        "exact examples include bra/ket coordinates under vectorisation, argmax/maxval of one parabola, " +
        "and the specific F81 Π-conjugation identity; inside/outside, classical/quantum and " +
        "inter-sectoral language remain interpretive comparisons. Operator-space vectorisation does not imply that " +
        "every object has exactly two readings, and the examples do not share polynomial ancestry.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode(
                "exact coordinate example: vectorisation",
                summary: "an operator on d dimensions has row and column indices and may be vectorised into d² coordinates; this fact alone does not classify all possible descriptions of the operator");
            yield return new InspectableNode(
                "exact scalar example: argmax and maxval",
                summary: "p(1-p) has argmax 1/2 and maxval 1/4; they are two different features of this named parabola, not a universal number/angle law");
            yield return new InspectableNode(
                "exact operator example: F81",
                summary: "F81 gives the specific identity Π·M·Π⁻¹ = M - 2L_{H_odd}; it relates two operator matrices without making every paired description an F81 descendant");
            yield return new InspectableNode(
                "interpretive catalogue",
                summary: "inside/outside calibration, classical/quantum reading modes and inter-sectoral reflection remain invitations to compare perspectives, not derived physical equivalences");
            yield return new InspectableNode(
                "open question",
                summary: "which additional pairs admit a precise common invariant, functor or change-of-coordinates statement that can fail under a counterexample?");
        }
    }
}
