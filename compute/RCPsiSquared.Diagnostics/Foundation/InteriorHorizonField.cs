using System.Globalization;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The Object Manager's telescope onto the third axis: the ¼-to-½ interior, read as a
/// horizon. Unlike the operator axes (crossover, J-defect), this is a state / coordinate axis: there
/// is no Hamiltonian to sweep, only the coherence CΨ approaching the cusp ¼. It makes the already
/// hardware-confirmed critical slowing (experiments/CRITICAL_SLOWING_AT_THE_CUSP.md) navigable live.
///
/// <para>Five readings, along a geometric |ε|-ladder approaching the horizon (CΨ = ¼ ± |ε|): the
/// marks (¼ the horizon, ½ the anchor); the heading θ → 0 from the interior (F95); the Mandelbrot
/// recursion run live, its iteration count diverging from the classical side (F56, the recursion
/// crawling at the fold, the horizon where time stops); the slowing-is-ours seam (a relative stop
/// makes the rescaled K constant, so the slowing belonged to the stop criterion, not the cusp); and
/// the γ-invariant dwell with the IBM Kingston anchors (F57). The horizon is structural (a saddle-node
/// fold), never gravitational.</para>
///
/// <para>The heading reads from the interior side (¼ + |ε|, plus the horizon point ¼ itself, where the
/// closed form is exactly θ = 0). The recursion reads from the classical side (¼ − |ε|, |ε| &gt; 0)
/// and never touches ¼: the Mandelbrot iteration has no fixed point at or beyond the cusp, so the
/// count is only defined below it.</para>
///
/// <para>A plain IInspectable, computed from closed forms plus the one live recursion. N-free (the
/// recursion and heading depend only on CΨ; state independence is itself hardware-confirmed).</para></summary>
public sealed class InteriorHorizonField : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private readonly double[] _eps;         // the |ε|-ladder, ascending
    private readonly double _tol;
    private readonly double _relK;
    private readonly double _gamma;

    public InteriorHorizonField(double epsLo = 1e-4, double epsHi = 0.25, int epsPoints = 13,
        double tol = 1e-12, double relK = 1e-3, double gamma = 0.5)
    {
        if (epsPoints < 2) throw new ArgumentOutOfRangeException(nameof(epsPoints), $"need at least two ε points; got {epsPoints}");
        if (epsLo <= 0 || epsHi <= epsLo || epsHi > 0.25)
            throw new ArgumentOutOfRangeException(nameof(epsHi), $"need 0 < epsLo < epsHi ≤ 0.25; got [{epsLo}, {epsHi}]");
        // These doubles flow straight into the closed forms (log/sqrt/divide), so the field owns their
        // validation: a non-positive tol, relK, or gamma would silently produce NaN or invert the grid.
        if (tol <= 0) throw new ArgumentOutOfRangeException(nameof(tol), $"tol must be positive; got {tol}");
        if (relK <= 0) throw new ArgumentOutOfRangeException(nameof(relK), $"relK must be positive; got {relK}");
        if (gamma <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), $"gamma must be positive; got {gamma}");
        _eps = GeometricLadder(epsLo, epsHi, epsPoints);
        _tol = tol;
        _relK = relK;
        _gamma = gamma;
    }

    private static double[] GeometricLadder(double lo, double hi, int count)
    {
        var grid = new double[count];
        double logLo = Math.Log(lo), logHi = Math.Log(hi);
        for (int i = 0; i < count; i++)
            grid[i] = Math.Exp(logLo + (logHi - logLo) * i / (count - 1));
        return grid;
    }

    /// <summary>The interior CΨ ladder for the heading (¼ + |ε|, ascending toward ½) with the horizon
    /// point ¼ prepended, where the heading is exactly θ = 0. The heading is a closed form valid at the
    /// cusp, so it may include the horizon itself; the recursion may not.</summary>
    private double[] HeadingCpsi()
    {
        var grid = new double[_eps.Length + 1];
        grid[0] = InteriorHorizon.Cusp;                          // the horizon, θ = 0 exactly
        for (int i = 0; i < _eps.Length; i++) grid[i + 1] = InteriorHorizon.Cusp + _eps[i];
        return grid;
    }

    /// <summary>The classical CΨ at each rung (¼ − |ε|), descending toward 0. Never includes ¼: the
    /// recursion has no fixed point at or beyond the cusp.</summary>
    private double[] ClassicalCpsi() => _eps.Select(e => InteriorHorizon.Cusp - e).ToArray();

    public string DisplayName => $"InteriorHorizonField (the ¼-to-½ interior, {_eps.Length} rungs, |ε| {_eps[0].ToString("E0", Inv)}..{_eps[^1].ToString("0.##", Inv)})";

    public string Summary
    {
        get
        {
            double thetaAnchor = InteriorHorizon.HeadingDegrees(InteriorHorizon.Cusp + _eps[^1]);
            int nNear = InteriorHorizon.RecursionIterations(InteriorHorizon.Cusp - _eps[0], _tol);
            return $"the horizon CΨ=¼: heading θ → 0 from the interior (θ={thetaAnchor.ToString("0.#", Inv)}° at the far rung), " +
                   $"the recursion crawls to {nNear} steps at the nearest rung (time stops at the fold); " +
                   "the slowing is ours (relative stop → constant); dwell γ-invariant (Kingston-confirmed). Structural horizon, no gravity.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // Readings are cheap closed forms plus one bounded live recursion (a few thousand steps
            // total across the rungs), so there is no per-enumeration cache here, unlike JDefectField
            // which caches because it rebuilds dense Liouvillians per sweep point.

            // 1. The marks (the contract): the horizon ¼ and the anchor ½, inert.
            yield return new InspectableNode(
                displayName: "the marks (the contract)",
                summary: $"CΨ=¼ the horizon (θ=0, the saddle-node cusp); CΨ=½ the anchor (θ=45°). Inert, tabulated, hardware-observed.");

            // 2. The heading from the interior: θ(CΨ) → 0 at the horizon (F95), the quantum side.
            var heading = HeadingCpsi();
            var theta = heading.Select(InteriorHorizon.HeadingDegrees).ToArray();
            yield return new InspectableNode(
                displayName: "the heading from the interior (θ → 0 at the horizon)",
                summary: $"θ = arctan(√(4·CΨ−1)): {theta[^1].ToString("0.#", Inv)}° at CΨ={heading[^1].ToString("0.##", Inv)} down to {theta[0].ToString("0.##", Inv)}° at the horizon. The compass on the quantum side (complex roots).",
                payload: new InspectablePayload.Curve("heading θ°", heading, theta, "CΨ (interior)", "θ°"));

            // 3. The recursion at the fold (live): the Mandelbrot iteration count diverges from below.
            var classical = ClassicalCpsi();
            var counts = classical.Select(c => (double)InteriorHorizon.RecursionIterations(c, _tol)).ToArray();
            double liveK = counts[0] * Math.Sqrt(_eps[0]);
            double closedK = InteriorHorizon.RecursionKClosedForm(classical[0], _tol);
            yield return new InspectableNode(
                displayName: "the recursion at the fold (the horizon, live)",
                summary: $"u → u²+c run live: {(int)counts[0]} steps at the nearest rung, diverging as CΨ → ¼⁻ (the recursion crawls, time stops). Rescaled K={liveK.ToString("0.##", Inv)} matches the closed form {closedK.ToString("0.##", Inv)}.",
                payload: new InspectablePayload.Curve("iteration count n", classical, counts, "CΨ (classical)", "n (live)"));

            // 4. The slowing is ours (the seam): relative stop → rescaled K constant.
            var kAbs = new double[_eps.Length];
            var kRel = new double[_eps.Length];
            for (int i = 0; i < _eps.Length; i++)
            {
                double c = classical[i];
                kAbs[i] = InteriorHorizon.RecursionIterations(c, _tol) * Math.Sqrt(_eps[i]);
                kRel[i] = InteriorHorizon.RecursionIterationsRelative(c, _relK) * Math.Sqrt(_eps[i]);
            }
            double kRelConst = 0.5 * Math.Log(4.0 / _relK);
            yield return new InspectableNode(
                displayName: "the slowing is ours (relative stop is flat)",
                summary: $"rescaled K with a relative stop tol=k·ε is the constant ½·ln(4/k)={kRelConst.ToString("0.##", Inv)}; the absolute-tol K drifts. The slowing belonged to the stop criterion, not the cusp (the cusp is inert).",
                children: new IInspectable[]
                {
                    new InspectableNode(
                        displayName: "K with absolute tol (drifts)",
                        summary: "the rescaled iteration count under a fixed tolerance",
                        payload: new InspectablePayload.Curve("K absolute", _eps, kAbs, "|ε|", "K")),
                    new InspectableNode(
                        displayName: "K with relative stop (flat = ½·ln(4/k))",
                        summary: $"converges to {kRelConst.ToString("0.###", Inv)}",
                        payload: new InspectablePayload.Curve("K relative", _eps, kRel, "|ε|", "K")),
                });

            // 5. The dwell and the hardware (the dose): the Bell+ geodesic + the Kingston anchors.
            var tGrid = GeometricTimeGrid();
            var geodesic = tGrid.Select(t => InteriorHorizon.BellPlusCpsi(_gamma, t)).ToArray();
            yield return new InspectableNode(
                displayName: "the dwell and the hardware (the dose)",
                summary: $"K_dwell = γ·t_dwell = {InteriorHorizon.BellPlusDwellPrefactor}·δ (F57), γ-invariant: a fixed dose carries Bell+ through the fold. " +
                         "Confirmed on IBM Kingston (f25_cusp_trajectory point-by-point; f57_kdwell_gamma_invariance two pairs, 6% spread, prefactor 0.67 under T1 damping). " +
                         "At the cusp θ→0 the Liouvillian eigenvalue is −γ₀ alone (pure decay): the horizon is where the carrier shows itself.",
                payload: new InspectablePayload.Curve("Bell+ geodesic CΨ(t)", tGrid, geodesic, "t", "CΨ (crosses ¼)"));
        }
    }

    private double[] GeometricTimeGrid()
    {
        // A time grid spanning the Bell+ crossing of the cusp (CΨ: 1/3 -> below 1/4). The crossing is
        // at K = gamma*t = 0.03735, i.e. t_cross = 0.03735/gamma; sample out to a few crossing times.
        double tCross = 0.03735 / _gamma;
        int points = 41;
        var grid = new double[points];
        for (int i = 0; i < points; i++) grid[i] = 4.0 * tCross * i / (points - 1);
        return grid;
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
