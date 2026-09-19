using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Coarse Q-regime band labels used by <see cref="QBasisAnker"/>. Mirrors the
/// three framework bands documented in <c>project_q_middle_structure</c> (onset, peak,
/// plateau) plus the Balance point (J = γ₀) and the candidate Endpoint-orbit position
/// outside the three core bands.</summary>
public enum QBand
{
    /// <summary>Onset band Q ∈ [0.2, 0.35]: dressed-mode-weight onset region.</summary>
    Onset,

    /// <summary>Balance point Q = 1: J = γ₀ exactly, the synchron of the two clocks
    /// (γ₀-Clock and H-Clock). See <c>hypotheses/Q_AS_THE_EXCHANGE_RATE.md</c>.</summary>
    Balance,

    /// <summary>Peak band Q ∈ [1.2, 1.8]: F86 Q_peak values per chromaticity sit here
    /// (c=2 → 1.5, c=3 → 1.6, c=4 and c=5 → 1.8 saturated).</summary>
    Peak,

    /// <summary>Plateau Q ≥ 2: idealised Q_EP at g_eff=1 sits here. Block-specific
    /// continuation past peak band.</summary>
    Plateau,

    /// <summary>Endpoint-orbit Q ≈ 2.5: stable across (c=2..4, N=5..8) per
    /// <see cref="F86.PerF71OrbitObservation"/>; candidate anchor not yet in
    /// QBasisAnkers proper.</summary>
    EndpointOrbit,
}

/// <summary>One typed Q-anchor on the Q = J/γ₀ axis. Carries the structural Q-value
/// plus four viewpoints: the corresponding J at the code-convention γ₀ = 0.05, the
/// band/role this anchor names, the tier label of its grounding, and the documenting
/// source.
///
/// <para><b>The 10 canonical anchors</b> (per <c>docs/Q_REGIME_ANCHORS.md</c>):
/// 0.2 onset start, 0.35 onset end, 1.0 Balance, 1.2 peak start, 1.5 F86 Q_peak c=2,
/// 1.6 F86 Q_peak c=3, sqrt(3) named two-level cross-reading, 1.8 F86 Q_peak c=4/c=5
/// + peak end, 2.0 Q_EP at g_eff=1 idealized, 2.5 Endpoint orbit (candidate).</para>
///
/// <para><b>Viewpoints</b> (Painter-Principle pluralism on Q):
/// <list type="bullet">
///   <item><see cref="Q"/>: the stored dimensionless coordinate for the named source.</item>
///   <item><see cref="JAtGamma0Point05"/>: the J = Q·γ₀ at the code convention.</item>
///   <item><see cref="Band"/>: which framework band this anchor names.</item>
///   <item><see cref="Role"/>: short prose description (e.g., "Balance", "F86 Q_peak (c=2)").</item>
///   <item><see cref="Tier"/>: typed knowledge tier per <see cref="Knowledge.Tier"/>.</item>
///   <item><see cref="DocumentingSource"/>: where the anchor is grounded.</item>
/// </list></para>
///
/// <para><b>Derived viewpoints</b> (helper methods): <see cref="JAt(double)"/> performs
/// the coordinate conversion J=Qγ₀; <see cref="ThetaDegrees"/> returns the generic dimensionless dial
/// arctan(Q) in degrees. A host may interpret that dial as an F95
/// principal angle only after it explicitly supplies the positive-decay roots
/// z=-lambda of a named lambda=-gamma0+/-iJ pair. It is not a universal many-body
/// mode angle.</para>
/// </summary>
public sealed record QBasisAnker(
    double Q,
    double JAtGamma0Point05,
    QBand Band,
    string Role,
    Tier Tier,
    string DocumentingSource)
{
    /// <summary>Coordinate conversion J = Q·γ₀ for the supplied γ₀.</summary>
    public double JAt(double gamma0) => Q * gamma0;

    /// <summary>The generic coordinate arctan(Q) in degrees. For the separately
    /// declared pair lambda=-gamma0+/-iJ with gamma0&gt;0 and J&gt;=0, its positive
    /// decay roots z=-lambda make this coordinate an F95 principal angle. Other
    /// Liouvillian modes need their own measured omega/gap. At Q=1 the coordinate
    /// is 45 degrees; as Q tends to infinity it tends to 90 degrees.</summary>
    public double ThetaDegrees() => AnchorConstants.RadiansToDegrees(Math.Atan(Q));
}
