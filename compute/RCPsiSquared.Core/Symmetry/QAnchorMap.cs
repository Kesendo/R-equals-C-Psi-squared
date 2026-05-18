using System.Text;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Typed Q-anchor map on the Q = J/γ₀ axis. The 9 canonical Q-anchors with
/// their viewpoints (Q, J at γ₀=0.05, band, role, tier, source), parallel to
/// <see cref="FractionReferenceGraph"/> on the α-axis and <see cref="PolarityMirrorMap"/>
/// on the γ-axis.
///
/// <para><b>Canonical 9 anchors</b> (per <c>docs/Q_REGIME_ANCHORS.md</c>):
/// onset edges 0.2 / 0.35, Balance 1.0, peak band 1.2 / 1.5 / 1.6 / 1.8, plateau
/// Q_EP 2.0, Endpoint orbit candidate 2.5.</para>
///
/// <para><b>What this map promotes</b>: the previous bare
/// <c>FractionReferenceGraph.QBasisAnkers = { 1.0, 1.5, 2.0 }</c> captured only the
/// three wave-breaking-scan anchors. The full 9-anchor structure was documented but
/// not typed; this map types it.</para>
///
/// <para><b>Wave-breaking-scan subset</b>: the original three QBasisAnkers (Q=1, 1.5,
/// 2.0) are the only anchors documented in <c>simulations/wave_breaking_q_anchor_scan.py</c>
/// with explicit dynamical-visibility properties (e.g., Q=1 hits F99 ankers at N=2,3,5).
/// The other six anchors (0.2, 0.35, 1.2, 1.6, 1.8, 2.5) are documented structural
/// markers without dedicated simulator sweeps.</para>
///
/// <para><b>What this map does NOT modify</b>: all F86 typed claims
/// (<see cref="F86.QEpLaw"/>, <see cref="F86.TPeakLaw"/>,
/// <see cref="F86.PerBlockQPeakClaim"/>, <see cref="F86.PerF71OrbitObservation"/>,
/// <c>BareDoubledPtfXPeak</c>, <c>HwhmLeftOverXPeakPrecise</c>, etc.) are referenced
/// but not touched. This map is an additive surface; the underlying F86 structure is
/// unchanged.</para>
///
/// <para><b>F86 integration</b>: surfaced as <c>F86KnowledgeBase.QAnchors</c> (block-
/// independent eager property, parallel to <c>PerBlockQPeaks</c>) and appears in the
/// <c>rcpsi inspect --root f86</c> tree under the top-level group "named Q-anchors
/// (Q = J/γ₀ axis)".</para>
/// </summary>
public sealed class QAnchorMap
{
    /// <summary>The 9 canonical Q-anchors. Constructed from the typed claims and
    /// documented Q-band edges they reference.</summary>
    public static IReadOnlyList<QBasisAnker> CanonicalAnchors { get; } =
        new QBasisAnker[]
        {
            new(Q: 0.2,
                JAtGamma0Point05: 0.010,
                Band: QBand.Onset,
                Role: "onset start",
                Tier: Tier.Tier2Empirical,
                DocumentingSource: "Q-band lower edge (project_q_middle_structure); ResonanceScan.DefaultQGrid lower bound"),

            new(Q: 0.35,
                JAtGamma0Point05: 0.0175,
                Band: QBand.Onset,
                Role: "onset end",
                Tier: Tier.Tier2Empirical,
                DocumentingSource: "Q-band upper edge (project_q_middle_structure)"),

            new(Q: 1.0,
                JAtGamma0Point05: 0.050,
                Band: QBand.Balance,
                Role: "Balance: J = γ₀ exactly",
                Tier: Tier.Tier1Derived,
                DocumentingSource: "wave_breaking_q_anchor_scan.py 'Balance (J=γ₀)'; synchron of γ₀-Clock and H-Clock"),

            new(Q: 1.2,
                JAtGamma0Point05: 0.060,
                Band: QBand.Peak,
                Role: "peak band start",
                Tier: Tier.Tier2Empirical,
                DocumentingSource: "Q-band lower edge (project_q_middle_structure)"),

            new(Q: 1.5,
                JAtGamma0Point05: 0.075,
                Band: QBand.Peak,
                Role: "F86 Q_peak (c=2); schema = Q_EP_central − r_polarity = 2 − 1/2 (Interior pole)",
                Tier: Tier.Tier1Derived,
                DocumentingSource: "PolarityPairQPeakDecompositionClaim (Tier1Derived schema 2 ± 1/2 inheriting QEpLaw + HalfAsStructuralFixedPoint); PerBlockQPeakClaim.Standard[0] empirical witness (wobble 1.4-1.6 across N=4..9, finite-size sensitive); wave_breaking_q_anchor_scan.py 'F86 Q_peak'"),

            new(Q: 1.6,
                JAtGamma0Point05: 0.080,
                Band: QBand.Peak,
                Role: "F86 Q_peak (c=3)",
                Tier: Tier.Tier2Empirical,
                DocumentingSource: "PerBlockQPeakClaim.Standard[1]; saturated N=5..9 (N-invariant)"),

            // The empirical F86 c=3 K-peak connection to √3 is open: fine-grid data
            // shows c=3 Interior Q_peak at N=7,8 = {1.743, 1.750} (within ~2% of √3),
            // but with non-monotone N-drift (+0.011 → +0.018), so the asymptote question
            // needs fine-grid scans at N=9,10,11+. This entry types only the bare-2×2
            // Lindblad-Absorption-Match anchor (Tier1Derived algebra).
            new(Q: 1.7320508075688772,
                JAtGamma0Point05: 0.05 * 1.7320508075688772,
                Band: QBand.Peak,
                Role: "Q=√3 canonical θ=60° Lindblad-Absorption-Match (|λ_±|=2γ₀=α)",
                Tier: Tier.Tier1Derived,
                DocumentingSource: "LindbladAbsorptionMatchAtSixtyDegreesClaim (Tier1Derived composition: F95 + AbsorptionTheorem + CanonicalTrigAnchor); bare-2×2 Lindblad eigenvalue magnitude equals single-site Absorption rate at the unique Q=√3, which is also a canonical Niven angle θ=60°; empirical F86 c=3 K-peak (N=7,8 ≈ 1.74-1.75) is suggestively close but the connection stays open (needs finer grid at N≥9)"),

            new(Q: 1.8,
                JAtGamma0Point05: 0.090,
                Band: QBand.Peak,
                Role: "F86 Q_peak (c=4, c=5) + peak band end",
                Tier: Tier.Tier2Empirical,
                DocumentingSource: "PerBlockQPeakClaim.Standard[2,3]; saturated N≥7; Q-band upper edge; the '1.8' is the coarse-grid Q_SCALE summary, fine-grid (PROOF_F86B_OBSTRUCTION) shows c=4 N=7=1.748, N=8=1.804 — structural asymptote stays open"),

            new(Q: 2.0,
                JAtGamma0Point05: 0.100,
                Band: QBand.Plateau,
                Role: "Q_EP at g_eff=1 (idealized)",
                Tier: Tier.Tier1Derived,
                DocumentingSource: "QEpLaw (Q_EP = 2/g_eff); wave_breaking_q_anchor_scan.py 'Q_EP (g_eff=1)'"),

            new(Q: 2.5,
                JAtGamma0Point05: 0.125,
                Band: QBand.EndpointOrbit,
                Role: "Endpoint orbit Q; schema = Q_EP_central + r_polarity = 2 + 1/2 (Endpoint pole)",
                Tier: Tier.Tier1Derived,
                DocumentingSource: "PolarityPairQPeakDecompositionClaim (Tier1Derived schema 2 ± 1/2 inheriting QEpLaw + HalfAsStructuralFixedPoint); PerF71OrbitObservation empirical witness (stable 2.39-2.61, ~2% N-variation across c=2..4 N=5..8); bit-exact value deviations Tier2Verified (PolarityInheritanceLink), closed-form blocked by PROOF_F86B_OBSTRUCTION"),
        };

    public IReadOnlyList<QBasisAnker> Anchors { get; }

    public QAnchorMap() : this(CanonicalAnchors) { }

    public QAnchorMap(IReadOnlyList<QBasisAnker> anchors)
    {
        Anchors = anchors ?? throw new ArgumentNullException(nameof(anchors));
    }

    /// <summary>The anchors in a specific band.</summary>
    public IReadOnlyList<QBasisAnker> ByBand(QBand band) =>
        Anchors.Where(a => a.Band == band).ToList();

    /// <summary>The anchors with a specific tier.</summary>
    public IReadOnlyList<QBasisAnker> ByTier(Tier tier) =>
        Anchors.Where(a => a.Tier == tier).ToList();

    /// <summary>Find the anchor with the given Q value (within tolerance). Returns
    /// null if no anchor matches.</summary>
    public QBasisAnker? AnkerAt(double q) =>
        Anchors.FirstOrDefault(a => Math.Abs(a.Q - q) < AnchorConstants.Tol);

    /// <summary>The three wave-breaking-scan anchors {1.0, 1.5, 2.0} — the subset
    /// with explicit dynamical-visibility properties documented in
    /// <c>simulations/wave_breaking_q_anchor_scan.py</c>.</summary>
    public IReadOnlyList<QBasisAnker> WaveBreakingScanSubset =>
        Anchors.Where(a => a.DocumentingSource.Contains("wave_breaking_q_anchor_scan")).ToList();

    /// <summary>Render the map as a printable table grouped by band.</summary>
    public string Render()
    {
        var sb = new StringBuilder();
        sb.AppendLine();
        sb.AppendLine("Q-Anchor Map (Q = J/γ₀ axis)");
        sb.AppendLine(new string('=', 88));
        sb.AppendLine();
        sb.AppendLine($"  Total anchors: {Anchors.Count}");
        sb.AppendLine($"  Wave-breaking-scan subset: {WaveBreakingScanSubset.Count}");
        sb.AppendLine();
        sb.AppendLine($"  {"Q",6} {"J(γ₀=0.05)",12} {"θ°",8} {"band",-16} {"tier",-22} {"role"}");
        sb.AppendLine($"  {new string('-', 6)} {new string('-', 12)} {new string('-', 8)} {new string('-', 16)} {new string('-', 22)} {new string('-', 30)}");

        foreach (var band in Enum.GetValues<QBand>())
        {
            var inBand = ByBand(band);
            if (inBand.Count == 0) continue;
            foreach (var a in inBand)
            {
                sb.AppendLine($"  {a.Q,6:F3} {a.JAtGamma0Point05,12:F4} {a.ThetaDegrees(),8:F1} {a.Band,-16} {a.Tier.Label(),-22} {a.Role}");
            }
        }

        sb.AppendLine();
        sb.AppendLine("  Wave-breaking-scan anchors (subset with documented dynamical visibility):");
        foreach (var a in WaveBreakingScanSubset)
        {
            sb.AppendLine($"    Q = {a.Q:F2}: {a.Role}");
        }

        return sb.ToString();
    }
}
