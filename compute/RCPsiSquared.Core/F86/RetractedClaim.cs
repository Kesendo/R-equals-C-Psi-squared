using RCPsiSquared.Core.Inspection;

using RCPsiSquared.Core.Knowledge;
namespace RCPsiSquared.Core.F86;

/// <summary>A formerly-claimed F86 closed form that was refuted by extended-N data. Kept
/// as a typed object — the OOP form of "intellectual honesty about the research process".
///
/// <para>Per the PTF retraction lesson (memory: <c>project_retraction_lesson_ptf_to_f86</c>):
/// closed-form conjectures from one (c, N) anchor can be trajectory crossings, not
/// asymptotes. What survives retraction is the symmetry, not the number. Encoding the
/// retracted claims here surfaces them in the OM tree so future investigations cannot
/// silently re-introduce them.</para>
/// </summary>
public sealed class RetractedClaim : Claim
{
    public string PreviousFormula { get; }
    public string Refutation { get; }

    public RetractedClaim(string name, string previousFormula, string refutation, string anchor)
        : base(name, Tier.Retracted, anchor)
    {
        PreviousFormula = previousFormula;
        Refutation = refutation;
    }

    public override string DisplayName => $"[RETRACTED] {Name}";
    public override string Summary => $"was: {PreviousFormula} → {Refutation}";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("previous formula", summary: PreviousFormula);
            yield return new InspectableNode("refutation", summary: Refutation);
        }
    }

    /// <summary>The retracted F86 closed forms / mechanisms: the two 2026-05-02 Q_peak
    /// csc(...) conjectures, plus the 2026-06-21 F86a "exceptional point on the real Q axis"
    /// mechanism (the Petermann factor is genuine non-normality on a simple eigenvalue, but its
    /// peak magnitude is a grid artifact). NOTE: the retraction's "no real-axis defective EP" was
    /// itself corrected after F89 certified narrow real-axis defective seeds at N=5,7,9. The
    /// all-odd endpoint-nullity surplus is not by itself an all-odd character theorem — see the
    /// Refutation string and PROOF_F86A_EP_MECHANISM §The real-axis EP.</summary>
    public static IReadOnlyList<RetractedClaim> Standard { get; } = new[]
    {
        new RetractedClaim(
            "Endpoint csc(π/(N+1))",
            "HWHM_left/Q_peak = csc(π/(N+1)) for Endpoint bonds",
            "N=7 coincidence; refuted by extended-N data (N=5..8 sweep diverges from formula)",
            "docs/proofs/PROOF_F86B_OBSTRUCTION.md (Retracted Q_peak conjectures, Endpoint) + memory: project_q_peak_ep_structure"),
        new RetractedClaim(
            "c=3 Interior csc(π/5)",
            "HWHM_left/Q_peak = csc(π/5) for c=3 Interior bonds",
            "N=7 coincidence; refuted by extended-N data",
            "docs/proofs/PROOF_F86B_OBSTRUCTION.md (Retracted Q_peak conjectures, c=3 Interior) + memory: project_q_peak_ep_structure"),
        new RetractedClaim(
            "F86a real-axis EP mechanism",
            "the (n,n+1)-coherence-block rate-channel degeneracy at Q_EP=2/g_eff is a defective exceptional point on the real Q axis; Petermann K ~6x FRAGILE_BRIDGE proves a real-axis EP",
            "The claimed Q_EP=2/g_eff mechanism is not a spectral degeneracy: the relevant eigenvalues there remain simple, while the large finite Petermann reading is genuine non-normality whose quoted peak magnitude was grid-sensitive. A separate, much narrower real-axis defective-seed family on the same block is certified by F89 at N=5,7,9 and was missed by the coarse F86a grid; the all-odd endpoint-nullity surplus does not alone extend that character verdict to every odd N. The Petermann peak law stays retracted; t_peak=1/(4gamma0) and Q_EP=2/g_eff only as the original channel-crossing definition survive.",
            "docs/proofs/PROOF_F86A_EP_MECHANISM.md (section The real-axis EP; verdict dates 2026-06-21 + 2026-07-07); typed at F86.LocalGlobalEpLink (OpenQuestion)"),
    };
}
