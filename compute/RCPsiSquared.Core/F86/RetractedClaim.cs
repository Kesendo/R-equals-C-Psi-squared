using RCPsiSquared.Core.Inspection;

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
public sealed class RetractedClaim : F86Claim
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

    /// <summary>The two retracted F86 closed forms from the 2026-05-02 retraction.</summary>
    public static IReadOnlyList<RetractedClaim> Standard { get; } = new[]
    {
        new RetractedClaim(
            "Endpoint csc(π/(N+1))",
            "HWHM_left/Q_peak = csc(π/(N+1)) for Endpoint bonds",
            "N=7 coincidence; refuted by extended-N data (N=5..8 sweep diverges from formula)",
            "docs/proofs/PROOF_F86_QPEAK.md retracted Statement 2 + memory: project_q_peak_ep_structure"),
        new RetractedClaim(
            "c=3 Interior csc(π/5)",
            "HWHM_left/Q_peak = csc(π/5) for c=3 Interior bonds",
            "N=7 coincidence; refuted by extended-N data",
            "docs/proofs/PROOF_F86_QPEAK.md retracted Statement 3 + memory: project_q_peak_ep_structure"),
    };
}
