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
    /// itself corrected 2026-07-07 (F89 locates a real-axis defective seed on this block at every
    /// odd N; the retraction's grid missed a √-EP window 20-30× narrower than its step) — see the
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
            "Independently re-computed (artifact-free) 2026-06-21: NO eigenvalue coalescence on the real Q axis - the relevant eigenvalues stay simple (nearest-neighbour gap ~0.25-0.35), so neither a defective EP nor a diabolic degeneracy sits there. The Petermann factor is GENUINE non-normality (the Riesz spectral-projector norm reproduces it on a simple isolated eigenvalue at Re~-4gamma0; N=5: ||P||=19.4=sqrt(375); cond(V)=49-268), NOT a degenerate-eigenspace eig artifact - but its peak MAGNITUDE is grid-sensitive (2-4x over deltaQ=1e-3), so '6x', '2384.7' and the within-parity growth law are grid artifacts, dropped. Adopted reading (PT_SYMMETRY_ANALYSIS): no real-axis EP; large-but-finite Petermann signals a nearby EP in the complex plane. Firmly-established defective EPs: the toy 2x2 reduction and the SEPARATE Sigma-gamma=0 gain-loss system (FRAGILE_BRIDGE). Whether the full Sigma-gamma=N*gamma0 block has an off-axis defective EP at all is OPEN (nearest complex-Q coalescences found are themselves diabolic). t_peak=1/(4gamma0) and Q_EP=2/g_eff-as-definition survive. SUPERSEDED IN PART 2026-07-07: F89's exact nullity count (r(0+)-r(inf)=N-1) proves a real-to-complex transition on this very (1,2) block at every odd N, and the Kato simple-zero lemma makes it defective (a Jordan block) at every seed tested, census-confirmed to N=11 and beta-exotic-scoped to Puiseux p~0.5 (only the codim-2 beta-exotic genericity stays open for all N). The 2026-06-21 real-axis scan missed it because a defective sqrt-EP splits its pair by ~sqrt|q-q*|, visible only within |q-q*|<~1e-3, while the scan's dQ~0.029 grid (121 pts over [0.5,4]) is 20-30x coarser and never sat inside that window (shown from below, F86aSeedMaskingTests). The Petermann-magnitude-as-grid-artifact call stands; the 'no real-axis defective EP' conclusion does not (superseded by F89's grid-robust count-change detector).",
            "docs/proofs/PROOF_F86A_EP_MECHANISM.md (section The real-axis EP; verdict dates 2026-06-21 + 2026-07-07); typed at F86.LocalGlobalEpLink (OpenQuestion)"),
    };
}
