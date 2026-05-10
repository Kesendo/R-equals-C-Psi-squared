namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>Thrown when a per-bond Q_peak in the c=2 K_b resonance scan lands within
/// one grid spacing of the upper Q-grid edge. This signals a methodological off-grid
/// escape (the true Q_peak lies beyond the scanned window) and means the reported
/// Q_peak/HWHM values for that bond are grid-snapped artefacts, not physical resonance
/// data.
///
/// <para>Tom Wicht 2026-05-10: "Lass uns auch eine Exception einbauen, wenn der Peak
/// exakt DefaultQGrid Grenze ist, das ist ja Fehlerquelle pur." The default grid
/// <see cref="Resonance.ResonanceScan.DefaultQGrid"/> is <c>LinearQGrid(0.20, 4.00, 153)</c>;
/// at N≥9 some flanking interior bonds peak above 4.0 and hit the edge silently. This
/// exception forces a conscious decision: extend the grid (pass a wider <c>qGrid</c> to
/// <see cref="C2HwhmRatio.Build"/>) or opt out of the validation (pass
/// <c>throwOnGridEdgeSnap: false</c>) when escape is the intended subject of study (see
/// <c>PerF71OrbitKTableTests.IsEscaped_FlagsFlankingOrbit_AtN9_WithDefaultGrid</c>).</para>
///
/// <para>Anchor: <c>memory/project_no_classicalization.md</c> — the off-grid escape is
/// methodological, NOT ontological: the orbit still has a Q_peak structure, we just
/// measured outside it. The exception keeps "we did not measure" distinct from "the
/// peak is at 4.0".</para></summary>
public sealed class GridEdgeEscapeException : InvalidOperationException
{
    /// <summary>The bond index whose Q_peak landed at the grid upper edge.</summary>
    public int Bond { get; }

    /// <summary>The grid-snapped Q_peak value (typically equal to the upper bound).</summary>
    public double QPeak { get; }

    /// <summary>The upper bound of the Q-grid that was used for the scan.</summary>
    public double GridUpper { get; }

    /// <summary>The grid spacing dQ at the upper edge (used for the within-one-spacing
    /// tolerance).</summary>
    public double DQ { get; }

    public GridEdgeEscapeException(int bond, double qPeak, double gridUpper, double dQ, string message)
        : base(message)
    {
        Bond = bond;
        QPeak = qPeak;
        GridUpper = gridUpper;
        DQ = dQ;
    }
}
