using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>F87 dissipator-axis selects polarity-axis (Tier 1 derived). Typed bridge between
/// <see cref="PolarityLayerOriginClaim"/> in <c>Pi2KnowledgeBase</c> (Core, "what is being
/// differentiated") and <see cref="DissipatorResonanceLaw"/> in <c>F87KnowledgeBase</c>
/// (Diagnostics, "where F87-hardness lives").
///
/// <para><b>Statement.</b> The dissipator letter is the choice of which polarity axis the
/// qubit's +0/−0 differentiation is read along. γ entering the cavity (cf.
/// <c>hypotheses/GAMMA_IS_LIGHT.md</c>, <c>hypotheses/RESONANCE_NOT_CHANNEL.md</c>) does
/// not merely decohere; it selects the axis on which the polarity becomes operational.
/// Z-dephasing activates the letter's own Klein cell (0, 1); X-dephasing activates
/// (1, 0); Y-dephasing activates (1, 1). The three letters are SU(2)-rotation-equivalent.</para>
///
/// <para><b>Not the Π² grading, and the distinction is load-bearing.</b> The cell above is
/// the KLEIN INDEX OF THE LETTER, which is where F87-hardness lives. It is a different map
/// from the parity Π_letter² grades by, which is bit_b for BOTH Z and Y and bit_a for X
/// alone (F38/F88a, typed at <see cref="PiOperator"/>, gated at
/// <c>simulations/f155_dephase_siblings.py</c> S5). Y's Klein index is (1, 1) while Π_Y²
/// grades by bit_b only; an earlier wording of this file read the two as one and asserted
/// that Π²_Y activates both axes at once, which is false. The conflation was NOT invented
/// here: the anchor <c>hypotheses/THE_POLARITY_LAYER.md</c> carried it in the same words
/// ("Z → bit_b, X → bit_a, Y → both", in a list explicitly of polarity AXES), and it is
/// repaired there in the same change.</para>
///
/// <para><b>Two readings unified.</b> The same multi-axis polarity layer
/// (Z₂² at k=2, Z₂³ at k≥3, per F88a) is read from <em>outside</em> by transverse-field
/// Brecher (h_y·Y or h_x·X breaking Z⊗N from the bit_a axis) or from <em>inside</em> by
/// dissipative resonance (F87-hardness in the matched Klein cell). Brecher and Hardness
/// are the two poles of dissipator-letter resonance; this Claim names the structural
/// fact that ties them.</para>
///
/// <para>Anchored in <c>hypotheses/THE_POLARITY_LAYER.md</c> §"Dissipator-resonance law"
/// and the 4×3 witness table at <see cref="DissipatorResonanceLaw.StandardWitnessTable"/>.
/// No new numbers are introduced; the bridge is structural, the witnesses live upstream.</para>
/// </summary>
public sealed class DissipatorAxisSelectsPolarityClaim : Claim
{
    public DissipatorAxisSelectsPolarityClaim()
        : base("dissipator letter selects polarity axis (bridge: PolarityLayerOrigin ↔ DissipatorResonance)",
               Tier.Tier1Derived,
               "hypotheses/THE_POLARITY_LAYER.md §Dissipator-resonance law + simulations/klein_dissipator_resonance.py + RCPsiSquared.Diagnostics/F87/DissipatorResonanceLaw + RCPsiSquared.Core/Symmetry/PolarityLayerOriginClaim")
    { }

    public override string DisplayName => "dissipator letter = polarity-axis selector (typed bridge)";

    public override string Summary =>
        "the dissipator letter selects which polarity axis the +0/−0 differentiation is read along; F87-hardness lives in the letter's own Klein cell, Z→(0,1), X→(1,0), Y→(1,1); that cell is NOT the Π² grading, which is bit_b for Z and Y and bit_a for X alone; SU(2)-equivalent; unifies Brecher (outside, transverse field) and Hardness (inside, F87) as two readings of one polarity";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("what is being differentiated",
                summary: "Pi2KnowledgeBase.PolarityLayerOrigin: multi-axis Z₂² (k=2) Klein-Vierergruppe with bit_a (X-eigenstate polarity) and bit_b (Π² parity); at k≥3 a third independent Y-parity axis (Z₂³, 8 sectors); +0/−0 lives on bit_a (cf. PolarityLayerOriginClaim layer 3)");
            yield return new InspectableNode("where F87-hardness lives",
                summary: "F87KnowledgeBase.DissipatorResonance: the 4×3 witness table: Mother (0,0) universally hard-free (0/66 across all letters); diagonal cells 50/76; off-diagonal 0/76; verified at N=4 k=3 (294 Z₂³-homogeneous pairs)");
            yield return new InspectableNode("the selector (Z → bit_b axis)",
                summary: "Z-dephasing has Klein index (0, 1); F87-hardness lives in Klein (0, 1); Π²_Z structure activates the bit_b axis of the polarity layer");
            yield return new InspectableNode("the selector (X → bit_a axis)",
                summary: "X-dephasing has Klein index (1, 0); F87-hardness lives in Klein (1, 0); Π²_X structure activates the bit_a axis (the +0/−0 axis itself becomes the dissipator-active axis)");
            yield return new InspectableNode("the selector (Y → Klein cell (1,1))",
                summary: "Y-dephasing has Klein index (1, 1); F87-hardness lives in Klein (1, 1). Its Π² grading is a SEPARATE map and is bit_b alone, the same axis as Π²_Z, never both at once (F38/F88a; gate S5 of simulations/f155_dephase_siblings.py). The polarity asymmetry follows that grading and not the Klein cell: it is exactly minus the Z value under Π_Y, since Π_Y = Π_Z⁻¹ (F155 §(g))");
            yield return new InspectableNode("two readings unified (Brecher ↔ Hardness)",
                summary: "Brecher = transverse field h_y·Y / h_x·X breaking Z⊗N from OUTSIDE the dissipator's Klein cell (bit_a-axis only); Hardness = F87 spectrum-pairing failure INSIDE the matched Klein cell, whose index depends on the letter and is (1,1) only for Y; same polarity layer, two perturbation types: unitary vs dissipative");
            yield return new InspectableNode("γ-as-light bridge",
                summary: "γ entering the cavity (hypotheses/GAMMA_IS_LIGHT.md) is not Shannon-noise; it is the choice of polarity-axis colour. The dissipator letter is the wavelength along which the +0/−0 differentiation becomes operational at d=2");
            yield return new InspectableNode("operational anchor",
                summary: "PauliPairTrichotomy.Classify(chain, terms, dephaseLetter: X|Y|Z) operationalises the axis-selection: the same Pauli pair classified under three different letters lands in three different Klein cells, hard in exactly one (the matched cell)");
        }
    }
}
