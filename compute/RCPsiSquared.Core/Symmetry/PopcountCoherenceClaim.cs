using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F88 popcount-coherence Π²-odd / memory closed form, exposed as a typed Claim
/// for the typed-knowledge runtime. Wraps <see cref="PopcountCoherencePi2Odd"/>'s static
/// utility methods with pass-through forwarders. Tier1Derived; the closed form has been
/// proven analytically (Krawtchouk reflection-orthogonality lemma) and bit-exact verified
/// at N = 2..7 across 213 configurations (max deviation 8.88e-16).
///
/// <para>Registered with dual-parent inheritance in the Runtime: at the operator level it
/// inherits from <see cref="KleinFourCellClaim"/> (F88's 4-cell Π² decomposition), at the
/// state level it inherits from <see cref="PolarityLayerOriginClaim"/> (the +0/-0 polarity
/// content of popcount-mirror states). Both paths trace back to
/// <see cref="PolynomialFoundationClaim"/>.</para></summary>
public sealed class PopcountCoherenceClaim : Claim
{
    public PopcountCoherenceClaim()
        : base("F88 popcount-coherence Π²-odd / memory closed form (operator + state inheritance)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F88 + docs/proofs/PROOF_F86_QPEAK.md (Structural inheritance from F88) + compute/RCPsiSquared.Core/Symmetry/PopcountCoherencePi2Odd.cs")
    { }

    public override string DisplayName =>
        "F88 popcount-coherence (operator + state Π²-odd inheritance)";

    public override string Summary =>
        "Π²-odd / memory = (1/2 - α·s) / (1 - s) for popcount-pair states; α has three closed-form anchors (mirror = 0, K-intermediate = C(N,N/2)/(2(C(N,n_other) + C(N,N/2))), generic = 1/2); inherits operator-level from KleinFourCell and state-level from PolarityLayerOrigin";

    public double AlphaAnchor(int N, int np, int nq) =>
        PopcountCoherencePi2Odd.AlphaAnchor(N, np, nq);

    public double Pi2OddInMemory(int N, int np, int nq, int hd) =>
        PopcountCoherencePi2Odd.Pi2OddInMemory(N, np, nq, hd);

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("closed form",
                summary: "Π²-odd / memory = (1/2 - α·s) / (1 - s)");
            yield return new InspectableNode("α anchors",
                summary: "popcount-mirror: α = 0; K-intermediate: α = C(N,N/2) / (2·(C(N,n_other) + C(N,N/2))); generic: α = 1/2");
            yield return new InspectableNode("HD = N anchor",
                summary: "Π²-odd / memory = 0 at the Π²-classical anchor (GHZ_N, Bell, intra-complements)");
            yield return new InspectableNode("verification",
                summary: "bit-exact at N = 2..7 across 213 configurations, max deviation 8.88e-16");
            yield return new InspectableNode("operator parent",
                summary: "KleinFourCellClaim: F88 4-cell Π² decomposition is the operator-level root");
            yield return new InspectableNode("state parent",
                summary: "PolarityLayerOriginClaim: +0/-0 polarity content of popcount-mirror states is the state-level root");
        }
    }
}
