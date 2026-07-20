using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;

using RCPsiSquared.Core.Knowledge;
namespace RCPsiSquared.Core.F86;

/// <summary>F86 Statement 1 partial: <c>Q_EP = 2/g_eff</c>, the definitional resonance Q,
/// where g_eff := σ_0, the SVD-top inter-channel singular value (no closed form). In the toy
/// 2-level reduction this is the genuine exceptional point at which the slowest eigenvalue pair
/// coalesces (PROOF_F86A_EP_MECHANISM Statement 1).
///
/// <para>Tier 1 derived in the 2-level reduction. <c>Q_EP = 2/g_eff</c> is a DEFINITION
/// (g_eff := σ_0); it is NOT "the full block-L exceptional-point position, bit-exact". The full
/// (n, n+1) block-L is genuinely NON-NORMAL near Q_peak but has NO defective EP on the real Q
/// axis; its eigenvalues stay simple there, and the large Petermann factor is finite and
/// grid-sensitive (the earlier "bit-exact full-block EP match" reading is retracted; see
/// <see cref="LocalGlobalEpLink"/>, now an OpenQuestion). For c=2, σ_0 → ≈ 2.8629 ± 1e-4 (the
/// F86e closure 2026-05-21, parity-split Aitken, γ-independent; SigmaZeroCommutatorNormClaim).
/// 2√2 was the N=7 finite-size crossing, not the limit (SigmaZeroChromaticityScaling demoted
/// it; doc line corrected 2026-06-10).</para>
/// </summary>
public sealed class QEpLaw : Claim
{
    public double GEff { get; }
    public double Value { get; }

    public QEpLaw(double gEff)
        : base("Q_EP law", Tier.Tier1Derived, "docs/ANALYTICAL_FORMULAS.md F86 Statement 1")
    {
        GEff = gEff;
        Value = EpAlgebra.QEp(gEff);
    }

    public override string DisplayName => "Q_EP = 2/g_eff";
    public override string Summary => $"= {Value:G6} (g_eff = {GEff:G6}, {Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("g_eff", GEff, "G6");
            yield return InspectableNode.RealScalar("Q_EP", Value, "G6");
        }
    }

    public override InspectablePayload Payload =>
        new InspectablePayload.Real("Q_EP", Value, "G6");
}
