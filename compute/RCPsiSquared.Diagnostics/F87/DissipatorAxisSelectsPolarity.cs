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
/// Z-dephasing activates Klein (0, 1) (bit_b axis); X-dephasing activates Klein (1, 0)
/// (bit_a axis); Y-dephasing activates Klein (1, 1) (both axes). The three letters are
/// SU(2)-rotation-equivalent.</para>
///
/// <para><b>Two readings unified.</b> The same multi-axis polarity layer
/// (Z₂² at k=2, Z₂³ at k≥3, per F88) is read from <em>outside</em> by transverse-field
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
        "the dissipator letter selects which polarity axis the +0/−0 differentiation is read along; Z→bit_b (Klein 0,1), X→bit_a (Klein 1,0), Y→both (Klein 1,1); SU(2)-equivalent; unifies Brecher (outside, transverse field) and Hardness (inside, F87) as two readings of one polarity";

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
            yield return new InspectableNode("the selector (Y → both axes)",
                summary: "Y-dephasing has Klein index (1, 1); F87-hardness lives in Klein (1, 1); Π²_Y activates bit_a + bit_b simultaneously");
            yield return new InspectableNode("two readings unified (Brecher ↔ Hardness)",
                summary: "Brecher = transverse field h_y·Y / h_x·X breaking Z⊗N from OUTSIDE the dissipator's Klein cell (bit_a-axis only); Hardness = F87 spectrum-pairing failure INSIDE the matched Klein cell (bit_a + bit_b); same polarity layer, two perturbation types: unitary vs dissipative");
            yield return new InspectableNode("γ-as-light bridge",
                summary: "γ entering the cavity (hypotheses/GAMMA_IS_LIGHT.md) is not Shannon-noise; it is the choice of polarity-axis colour. The dissipator letter is the wavelength along which the +0/−0 differentiation becomes operational at d=2");
            yield return new InspectableNode("operational anchor",
                summary: "PauliPairTrichotomy.Classify(chain, terms, dephaseLetter: X|Y|Z) operationalises the axis-selection: the same Pauli pair classified under three different letters lands in three different Klein cells, hard in exactly one (the matched cell)");
        }
    }
}
