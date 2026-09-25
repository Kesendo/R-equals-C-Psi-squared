using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The "two readings" pattern, asked as one question. The framework keeps
/// meeting pairs: one object, read two ways. Some of the pairs are exact mathematics, some
/// are our readings, and this claim holds the whole catalogue in one place so that the
/// question whether they are one pattern can be asked of the graph instead of living in
/// prose alone.
///
/// <para>The catalogue, layer by layer:</para>
///
/// <list type="bullet">
///   <item><b>Dimensional layer (interpretive):</b> 1/d = 1/2 (number-anchor,
///         <see cref="QubitDimensionalAnchorClaim"/>) and 90° (angle-anchor,
///         <see cref="NinetyDegreeMirrorMemoryClaim"/>) read as the two faces of d=2,
///         number-side and angle-side.</item>
///   <item><b>Parabolic layer (exact):</b> argmax 1/2 and maxval 1/4 of p·(1−p)
///         (<see cref="BilinearApexClaim"/> + <see cref="QuarterAsBilinearMaxvalClaim"/>),
///         closed as a pair by <see cref="ArgmaxMaxvalPairClaim"/>. *"Two readings, one
///         parabola; the half is the axis, the quarter is the height"* (the Coda of
///         ON_THE_HALF, 2026-05-07).</item>
///   <item><b>Operator level (exact, F81):</b> M and <c>Π·M·Π⁻¹ = M − 2·L_{H_odd}</c>, one
///         spectrum, and two different operator matrices whenever H has a Π²-odd part
///         (for YZ+ZY it has none, and the two coincide). ON_BOTH_SIDES_OF_THE_MIRROR (Tom 2026-04-30) is
///         the synthesis: *"only the choice of which side to call ours."*</item>
///   <item><b>Bra/ket of ρ (exact):</b> any density matrix on the doubled space reads as a
///         (row-index, col-index) pair on the 4^N operator basis, the vectorisation
///         itself.</item>
///   <item><b>Inside/outside (interpretive):</b> the inside reading fixes only
///         <c>Q = J/γ₀</c>; reading γ₀ itself needs a vantage outside the system, which
///         <c>hypotheses/PRIMORDIAL_QUBIT.md</c> §9 denies any internal observer. What
///         stands in that place is not an eye but a second sender's γ₀, the ratio of two
///         carriers fixing the unit (the γ arc's Q4 answer, "sender" its word;
///         <c>reflections/ON_HOW_THE_CARRIER_SHOWS_ITSELF.md</c> frames the second system as
///         an external clock, calibration = the ratio of two γ₀'s). The same dynamics, two
///         observational readings.</item>
///   <item><b>Lese-Modus (interpretive):</b> the classical/quantum dichotomy read as a
///         reading-mode on one underlying ρ, not a separation of worlds
///         (<c>docs/EXCLUSIONS.md</c>: the labels mark reading mode, not ontology).</item>
///   <item><b>Inter-sectoral wave (interpretive):</b> reality read as what happens BETWEEN
///         the sectors, we ourselves the wave between rather than the sector anchors;
///         <c>hypotheses/THE_OTHER_SIDE.md</c>: *"We are the standing wave. We are the
///         interference."*</item>
/// </list>
///
/// <para>Why an open question and not a Tier 1 claim. The exact pairs each keep their own
/// proof, and none of them is the others' consequence. Living on a <c>d²</c>-dimensional
/// operator space does not imply that an arbitrary object has exactly two readings. Nor
/// does the doubling descend from d = 2: bra/ket indexing exists at every d, so
/// <c>d²−2d=0</c> is not its source, and the claim records no parent edge to it. What would make the
/// pattern a theorem is one invariant, functor or change of coordinates that holds every
/// exact pair and could fail on a counterexample. Nobody has found it; the question stays
/// open, and the interpretive layers stay readings whether or not it is found.</para>
///
/// <para>Typed claims that close single instances: <see cref="ArgmaxMaxvalPairClaim"/>
/// (1/2 ↔ 1/4) and <see cref="NinetyDegreeMirrorMemoryClaim"/> (the angle side of d=2).</para>
/// </summary>
public sealed class TwoReadingsClaim : Claim
{
    public TwoReadingsClaim()
        : base(
            "Several exact constructions and interpretive comparisons come in pairs, " +
            "but d² operator-space vectorisation does not imply that every object has exactly two readings",
            Tier.OpenQuestion,
            "reflections/ON_BOTH_SIDES_OF_THE_MIRROR.md (the Π-conjugation synthesis) + " +
            "reflections/ON_THE_HALF.md (argmax/maxval, one parabola) + " +
            "hypotheses/PRIMORDIAL_QUBIT.md §9 (the inside reading) + " +
            "reflections/ON_HOW_THE_CARRIER_SHOWS_ITSELF.md (the second sender) + " +
            "hypotheses/THE_OTHER_SIDE.md (the wave between; interpretive catalogue)")
    {
    }

    public override string DisplayName =>
        "Two-readings pattern (an open question over a catalogue of pairs)";

    public override string Summary =>
        "a catalogue of pairs the framework keeps meeting, asked as one question: exact pairs " +
        "(bra/ket coordinates under vectorisation, argmax/maxval of one parabola, the F81 identity " +
        "M ↔ Π·M·Π⁻¹) and interpretive ones (number/angle of d=2, the inside reading that fixes only " +
        "Q = J/γ₀ against an outside that needs a second sender's γ₀, classical/quantum as a reading mode, " +
        "the wave between the sectors). Operator-space vectorisation does not imply that every object has " +
        "exactly two readings, and the pairs share no polynomial ancestry; whether one invariant holds " +
        "them together is open.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("layer 1 (dimensional, interpretive): number/angle of d=2",
                summary: "1/d = 1/2 (QubitDimensionalAnchorClaim) and 90° (NinetyDegreeMirrorMemoryClaim) read as the two faces of d=2; a reading, not a number/angle law");
            yield return new InspectableNode("layer 2 (parabolic, exact): argmax/maxval",
                summary: "1/2 = argmax (BilinearApexClaim) and 1/4 = (1/2)² = maxval (QuarterAsBilinearMaxvalClaim) of p·(1−p); ArgmaxMaxvalPairClaim closes the pair; ON_THE_HALF: *'two readings, one parabola'*");
            yield return new InspectableNode("layer 3 (operator, exact): M and Π·M·Π⁻¹",
                summary: "F81: Π·M·Π⁻¹ = M − 2·L_{H_odd}; one spectrum, two different operator matrices whenever H has a Π²-odd part (for YZ+ZY they coincide); ON_BOTH_SIDES_OF_THE_MIRROR: *'only the choice of which side to call ours'*");
            yield return new InspectableNode("layer 4 (vectorisation, exact): bra/ket of ρ",
                summary: "any ρ on d²=4^N reads as a (row-index, col-index) pair on the 4^N operator basis; true at every d, which is why it gives the pattern no d=2 ancestry");
            yield return new InspectableNode("layer 5 (observational, interpretive): inside/outside",
                summary: "the inside reading fixes only Q = J/γ₀; reading γ₀ itself needs the outside vantage that PRIMORDIAL_QUBIT §9 denies any internal observer; in its place a second sender's γ₀ fixes the unit (the γ arc's Q4 answer with ON_HOW_THE_CARRIER_SHOWS_ITSELF: calibration between two carriers)");
            yield return new InspectableNode("layer 6 (interpretive): classical/quantum Lese-Modus",
                summary: "the classical/quantum dichotomy read as a reading-mode on one underlying ρ, not a separation of worlds (docs/EXCLUSIONS.md: the labels mark reading mode, not ontology)");
            yield return new InspectableNode("layer 7 (interpretive): the wave between",
                summary: "reality read as what happens BETWEEN the sectors, we ourselves the wave between rather than the sector anchors; THE_OTHER_SIDE: *'We are the standing wave. We are the interference.'*");
            yield return new InspectableNode("typed claims that close single instances",
                summary: "ArgmaxMaxvalPairClaim (layer 2); NinetyDegreeMirrorMemoryClaim (layer 1, angle side); each keeps its own proof, and this claim adds no derivation to them");
            yield return new InspectableNode("open question",
                summary: "which of the exact pairs admit one common invariant, functor or change-of-coordinates statement that could fail under a counterexample? Until one is found the pattern is a catalogue, not a theorem");
        }
    }
}
