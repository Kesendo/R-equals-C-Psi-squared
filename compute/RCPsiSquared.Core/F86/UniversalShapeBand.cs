using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86;

/// <summary>The F86 Statement 2 band per bond class: the smallest and the largest envelope
/// HWHM_left/Q_peak among the class's witnesses that have been read at the F86 standard locations
/// (<see cref="F86StandardLocations"/>), by <see cref="ResonanceScan"/> at γ₀ = 0.05 on the default
/// Q grid. The prediction's expected ratio is the band's midpoint and its tolerance the half-width
/// plus <see cref="Resolution"/>, so every read witness lies inside by construction: the band is a
/// record of the span, not a prediction that could fail. Every standard location is read; the
/// Interior values rise with N at each c (c=3: 0.7572 → 0.7622 over N=5..8; c=4: 0.7623 → 0.7637
/// over N=7, 8), so a witness at larger N may lie above <see cref="InteriorHigh"/>.
///
/// <para>Second route for the edges InteriorLow, EndpointLow and EndpointHigh:
/// simulations/f86_envelope_second_route.py (the fixture read by
/// ValidateAgainstPythonEnvelopeRouteTests), which reproduces the grid reading to 10⁻¹⁰.
/// InteriorHigh (c=4 N=8, block dimension 3920) is read by ResonanceScan alone so far; the
/// script reads it with `--fixture 8 4`.
/// The grid reading itself is off the continuous envelope peak by up to 3.5·10⁻⁵ in the ratio per
/// bond at c=2 N=5..8 (the class-mean edges were not checked against a continuous search), from the parabolic vertex and the linear half-maximum on a grid of step 0.025
/// (simulations/f86_envelope_second_route.py --continuous N), so on another Q grid the edges move
/// by that much.</para></summary>
public static class UniversalShapeBand
{
    /// <summary>Interior low edge: c=2 N=5.</summary>
    public const double InteriorLow = 0.7523433331317626;

    /// <summary>Interior high edge: c=4 N=8.</summary>
    public const double InteriorHigh = 0.7636747635673543;

    /// <summary>Endpoint low edge: c=3 N=5.</summary>
    public const double EndpointLow = 0.7710090593327815;

    /// <summary>Endpoint high edge: c=2 N=6.</summary>
    public const double EndpointHigh = 0.7815203313096158;

    /// <summary>The scan's reproducibility in the ratio on one grid: a live re-read of an edge
    /// witness differs from the frozen literal by a few 10⁻¹¹ (the peak value K is fixed to ~10⁻¹³
    /// by the eigensolver's round-off, and the parabolic vertex on a Q grid of step 0.025 magnifies
    /// it by 1/(|curvature|·step)), so the band is widened by this much on each side. It is not the
    /// reading's accuracy (see the class summary).</summary>
    public const double Resolution = 1e-10;

    public static UniversalShapePrediction Interior(IReadOnlyList<UniversalShapeWitness> witnesses) =>
        Build(BondClass.Interior, InteriorLow, InteriorHigh, witnesses);

    public static UniversalShapePrediction Endpoint(IReadOnlyList<UniversalShapeWitness> witnesses) =>
        Build(BondClass.Endpoint, EndpointLow, EndpointHigh, witnesses);

    private static UniversalShapePrediction Build(BondClass cls, double low, double high,
        IReadOnlyList<UniversalShapeWitness> witnesses) =>
        new(cls, expectedRatio: 0.5 * (low + high), tolerance: 0.5 * (high - low) + Resolution, witnesses);
}
