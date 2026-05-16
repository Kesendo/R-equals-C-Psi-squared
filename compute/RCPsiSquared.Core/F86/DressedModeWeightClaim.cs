using RCPsiSquared.Core.Inspection;

using RCPsiSquared.Core.Knowledge;
namespace RCPsiSquared.Core.F86;

/// <summary>F86 structural fingerprint: at Q = Q_peak the Dicke probe sits dominantly in
/// dressed (H-mixed, complex-eigenvalue) modes; far past Q_peak it sits substantially in
/// pure-rate modes again. This is what makes Q_peak a generalised exceptional-point
/// resonance condition: ∂S/∂J peaks where probe weight has been pulled off the pure-rate
/// ladder onto the first complex-conjugate eigenvalue pair just past the EP.
///
/// <para><b>Tier 1 candidate.</b> The qualitative structural mechanism (probe-weight
/// shift onto dressed pair past EP) is sound — it follows from the F86a EP analysis.
/// The specific hardcoded values <see cref="WeightAtQPeak"/> = 0.99 and
/// <see cref="WeightAtPlateau"/> = 0.31 at <see cref="PlateauQ"/> = 20 are
/// <b>unverified anchors</b>, not derived constants and not witness-backed per (c, N):</para>
///
/// <list type="bullet">
///   <item>The empirical W_peak table in <c>experiments/Q_SCALE_THREE_BANDS.md</c>
///         (lines 95–123) shows W_peak ranges 0.832 (N=4 c=2) → 0.9996 (N=9 c=3); 0.99
///         is the high-c high-N asymptote, NOT a universal constant.</item>
///   <item>The empirical W_plateau is measured at Q = 50 (not Q = 20) in the same table
///         and ranges 0.42–0.86 across blocks; the hardcoded 0.31 sits below that
///         empirical range, indicating the value comes from a DIFFERENT Q point or a
///         different (c, N) than the canonical table.</item>
///   <item>No C# code computes W(Q) per (c, N); this Claim is two memorised constants
///         with no per-block witness collection.</item>
/// </list>
///
/// <para><b>To promote Tier 1 candidate → Tier 1 derived:</b> either (a) derive
/// W(Q_peak) and W(Q_plateau) analytically from the EP 2-level eigenvector rotation, or
/// (b) replace the hardcoded constants with a per-(c, N) witness collection backed by
/// actual computation (`simulations/_eq022_b1_step_c_time_evolution.py` has the
/// `dressed_weight` calculator). Per-block tracking would also distinguish c-dependence
/// from N-dependence.</para></summary>
public sealed class DressedModeWeightClaim : Claim
{
    /// <summary>Approximate W at Q_peak from one historical simulation; high-c high-N
    /// asymptote per <c>experiments/Q_SCALE_THREE_BANDS.md</c>, NOT a universal constant
    /// (W_peak ranges 0.832 at N=4 c=2 to 0.9996 at N=9 c=3).</summary>
    public double WeightAtQPeak { get; } = 0.99;

    /// <summary>Hardcoded W at <see cref="PlateauQ"/>; sits below the empirical W_plateau
    /// range (0.42–0.86 at Q=50 per Q_SCALE_THREE_BANDS.md). Likely from a different Q
    /// point or a specific small-N c=2 case; not witness-backed per (c, N).</summary>
    public double WeightAtPlateau { get; } = 0.31;

    /// <summary>Q at which <see cref="WeightAtPlateau"/> was measured; canonical
    /// Q_SCALE_THREE_BANDS plateau is at Q = 50, NOT 20.</summary>
    public double PlateauQ { get; } = 20.0;

    public DressedModeWeightClaim()
        : base("dressed-mode probe weight at Q_peak (qualitative structural fingerprint; specific values unverified anchors)",
               Tier.Tier1Candidate,
               "docs/ANALYTICAL_FORMULAS.md F86 EP-based structural mechanism (qualitative); " +
               "experiments/Q_SCALE_THREE_BANDS.md (empirical W tables, NOT matching hardcoded 0.99/0.31); " +
               "simulations/_eq022_b1_step_c_time_evolution.py (dressed_weight calculator)")
    { }

    public override string DisplayName =>
        "Dicke probe → dressed modes at Q_peak (Tier 1 candidate; structural mechanism, specific values unverified)";

    public override string Summary =>
        $"qualitative: W(Q_peak) >> W(Q_plateau); hardcoded anchors W(Q_peak)≈{WeightAtQPeak:P0}, W(Q={PlateauQ:F0})≈{WeightAtPlateau:P0} unverified per (c, N) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("dressed-mode weight at Q_peak (anchor, unverified)", WeightAtQPeak, "P1");
            yield return InspectableNode.RealScalar($"dressed-mode weight at Q={PlateauQ} (anchor, unverified)", WeightAtPlateau, "P1");
            yield return new InspectableNode("structural mechanism (Tier 1)",
                summary: "probe pulled off pure-rate ladder onto complex-conjugate pair just past EP — Q_peak is a generalised EP resonance");
            yield return new InspectableNode("specific values caveat",
                summary: "0.99/0.31 are memorised anchors; empirical W_peak in Q_SCALE_THREE_BANDS.md ranges 0.832-0.9996, W_plateau at Q=50 ranges 0.42-0.86. No per-(c, N) witness collection.");
        }
    }
}
