using System.Globalization;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The Object Manager's telescope onto the third axis: the ¼-to-½ interior, read as a
/// horizon. Unlike the operator axes (crossover, J-defect), this is a coordinate axis: there is no
/// Hamiltonian to sweep, only the recurrence value c approaching its double root ¼. It shows the
/// recurrence discriminant, the F95 heading at b = ½, the live iteration count, its stopping rule,
/// and one named Bell+/pure-Z trajectory (F25) that passes through the same scalar value. The
/// horizon is structural, never gravitational.
///
/// <para>Five readings use a geometric |ε|-ladder around ¼: the scalar marks; the F95 heading
/// θ → 0 from c &gt; ¼; the recurrence iteration count from c &lt; ¼; the relative-stop control
/// showing how much of the rescaled count belongs to the stopping rule; and the Bell+/pure-Z dwell
/// with its Kingston records (experiments/CRITICAL_SLOWING_AT_THE_CUSP.md).</para>
///
/// <para>The heading reads from the complex-root side (¼ + |ε|, plus ¼ itself where the
/// closed form is exactly θ = 0). The recurrence reads from the two-real-root side
/// (¼ − |ε|, |ε| &gt; 0) and never touches the double-root point.</para>
///
/// <para>A plain IInspectable, computed from closed forms, the F56 asymptotic, and one live recursion. N-free (the
/// recurrence and angle depend only on their scalar input; that is not a state-independence claim).</para></summary>
public sealed class InteriorHorizonField : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The near end of the default |ε|-ladder. Named rather than repeated, so the CLI's
    /// <c>--eps-lo</c> fallback is this value and not a second copy of it.</summary>
    public const double DefaultEpsLo = 1e-10;

    /// <summary>The far end of the default |ε|-ladder, kept inside F56's committed scan range.</summary>
    public const double DefaultEpsHi = 1e-2;

    private readonly double[] _eps;         // the |ε|-ladder, ascending
    private readonly double _tol;
    private readonly double _relK;
    private readonly double _gamma;

    /// <summary>The default ladder runs 10⁻¹⁰ … 10⁻² because that is inside F56's committed scan
    /// range (tol 10⁻⁸…10⁻¹⁶, ε 10⁻¹…10⁻¹⁰) and inside its stated validity <c>tol ≪ ε ≪ 1</c>. The
    /// constructor still accepts ε up to ¼, but a rung there is not a useful near-boundary reading: at
    /// CΨ = ¼ − ¼ = 0 the iteration starts at u₀ = 0, its first increment is 0, and the count is 1
    /// whatever the physics does.</summary>
    public InteriorHorizonField(double epsLo = DefaultEpsLo, double epsHi = DefaultEpsHi, int epsPoints = 13,
        double tol = 1e-12, double relK = 1e-3, double gamma = 0.5)
    {
        if (epsPoints < 2) throw new ArgumentOutOfRangeException(nameof(epsPoints), $"need at least two ε points; got {epsPoints}");
        if (epsLo <= 0) throw new ArgumentOutOfRangeException(nameof(epsLo), $"epsLo must be positive; got {epsLo}");
        if (epsHi <= epsLo || epsHi > 0.25)
            throw new ArgumentOutOfRangeException(nameof(epsHi), $"need epsLo < epsHi ≤ 0.25; got [{epsLo}, {epsHi}]");
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

    /// <summary>The c&gt;¼ ladder for the angle (¼ + |ε|, ascending toward ½), with the
    /// double-root point ¼ prepended where the angle is exactly zero.</summary>
    private double[] HeadingCpsi()
    {
        var grid = new double[_eps.Length + 1];
        grid[0] = InteriorHorizon.Cusp;                          // recurrence boundary, θ = 0 exactly
        for (int i = 0; i < _eps.Length; i++) grid[i + 1] = InteriorHorizon.Cusp + _eps[i];
        return grid;
    }

    /// <summary>The two-real-root recurrence inputs at each rung (¼ − |ε|), descending toward 0.
    /// Never includes the double-root point ¼.</summary>
    private double[] BelowBoundaryCpsi() => _eps.Select(e => InteriorHorizon.Cusp - e).ToArray();

    public string DisplayName => $"InteriorHorizonField (recurrence coordinate near ¼, {_eps.Length} rungs, |ε| {_eps[0].ToString("E0", Inv)}..{_eps[^1].ToString("0.##", Inv)})";

    public string Summary
    {
        get
        {
            double thetaAnchor = InteriorHorizon.HeadingDegrees(InteriorHorizon.Cusp + _eps[^1]);
            int nNear = InteriorHorizon.RecursionIterations(InteriorHorizon.Cusp - _eps[0], _tol);
            return $"recurrence boundary c=¼: the F95 heading (b=½) θ → 0 from the complex-root side (θ={thetaAnchor.ToString("0.#", Inv)}° at the far rung), " +
                   $"the recurrence takes {nNear} steps at the nearest two-real-root rung under the fixed stop; " +
                   "a relative stop makes the rescaled count approach a constant, while the Bell+/pure-Z dwell is a separate setup-specific reading. A structural horizon, not a gravitational one.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // Readings are cheap closed forms plus one bounded live recursion (a few thousand steps
            // total across the rungs), so there is no per-enumeration cache here, unlike JDefectField
            // which caches because it rebuilds dense Liouvillians per sweep point.

            // 1. The scalar marks: the recurrence boundary ¼ and the angle anchor ½.
            yield return new InspectableNode(
                displayName: "the marks (the contract)",
                summary: $"c=¼ is the recurrence's double-root value and θ=0; c=½ is an angle anchor with θ=45°. These are scalar coordinates, not state-regime labels.");

            // 2. The F95 heading at b = ½ on the complex-root side.
            var heading = HeadingCpsi();
            var theta = heading.Select(InteriorHorizon.HeadingDegrees).ToArray();
            yield return new InspectableNode(
                displayName: "the F95 heading at b = ½ (θ → 0 at the double root)",
                summary: $"θ = arctan(√(4·c−1)) = F95's θ(c; ½) of z²−z+c: {theta[^1].ToString("0.#", Inv)}° at c={heading[^1].ToString("0.##", Inv)} down to {theta[0].ToString("0.##", Inv)}° at the double root ¼. The angle of the recurrence's complex root pair.",
                payload: new InspectablePayload.Curve("heading θ°", heading, theta, "c (complex-root side)", "θ°"));

            // 3. The recurrence at its double-root boundary, run live from below.
            var belowBoundary = BelowBoundaryCpsi();
            var counts = belowBoundary.Select(c => (double)InteriorHorizon.RecursionIterations(c, _tol)).ToArray();
            double liveK = counts[0] * Math.Sqrt(_eps[0]);
            double asymptoticK = InteriorHorizon.RecursionKAsymptotic(belowBoundary[0], _tol);
            yield return new InspectableNode(
                displayName: "the recurrence at the double-root boundary (live)",
                summary: $"u → u²+c run live: {(int)counts[0]} steps at the nearest rung under the fixed stopping rule as c → ¼⁻. Rescaled K={liveK.ToString("0.##", Inv)}; the finite-ε asymptotic gives {asymptoticK.ToString("0.##", Inv)}.",
                payload: new InspectablePayload.Curve("iteration count n", belowBoundary, counts, "c (two-real-root side)", "n (live)"));

            // 4. The slowing is ours (the seam): relative stop → rescaled K constant.
            var kAbs = new double[_eps.Length];
            var kRel = new double[_eps.Length];
            for (int i = 0; i < _eps.Length; i++)
            {
                double c = belowBoundary[i];
                kAbs[i] = InteriorHorizon.RecursionIterations(c, _tol) * Math.Sqrt(_eps[i]);
                kRel[i] = InteriorHorizon.RecursionIterationsRelative(c, _relK) * Math.Sqrt(_eps[i]);
            }
            double kRelConst = 0.5 * Math.Log(4.0 / _relK);
            yield return new InspectableNode(
                displayName: "the slowing is ours (relative stop flattens onto ½·ln(4/k))",
                summary: $"rescaled K with a relative stop tol=k·ε tends to ½·ln(4/k)={kRelConst.ToString("0.###", Inv)} as ε → 0, reading {kRel[0].ToString("0.###", Inv)} at the nearest rung |ε|={_eps[0].ToString("E0", Inv)} and drifting to {kRel[^1].ToString("0.###", Inv)} at |ε|={_eps[^1].ToString("E0", Inv)}, where ε is no longer small; the absolute-tol K drifts across the whole ladder ({kAbs[0].ToString("0.##", Inv)} to {kAbs[^1].ToString("0.##", Inv)}). This control separates the stopping-rule contribution from the recurrence boundary.",
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

            // 5. A separate Bell+/pure-Z trajectory and finite Kingston comparisons.
            var tGrid = GeometricTimeGrid();
            var trajectory = tGrid.Select(t => InteriorHorizon.BellPlusCpsi(_gamma, t)).ToArray();
            yield return new InspectableNode(
                displayName: "the named Bell+/pure-Z dwell and hardware comparison",
                summary: $"K_dwell = γ·t_dwell = {InteriorHorizon.BellPlusDwellPrefactor}·δ (F57) for the ideal Bell+/pure-Z trajectory. " +
                         "On IBM Kingston, f25_cusp_trajectory reproduces the F25 crossing point by point (19 delay points, RMS residual 0.0097 against the in-situ γ fit); " +
                         "f57_kdwell_gamma_invariance compares two pairs at 2.55× different fitted γ, their 6.4% agreement at prefactor 0.67 rather than 1.0801 not isolating the cause of the gap. " +
                         "Neither the scalar crossing nor that association identifies a recurrence root with a Liouvillian mode.",
                payload: new InspectablePayload.Curve("Bell+/pure-Z scalar trajectory", tGrid, trajectory, "t", "CΨ (passes ¼)"));
        }
    }

    private double[] GeometricTimeGrid()
    {
        // A time grid spanning the Bell+ scalar crossing (CΨ: 1/3 -> below 1/4). The crossing is
        // at K = gamma*t = 0.03735, i.e. t_cross = 0.03735/gamma; sample out to a few crossing times.
        double tCross = 0.03735 / _gamma;
        int points = 41;
        var grid = new double[points];
        for (int i = 0; i < points; i++) grid[i] = 4.0 * tCross * i / (points - 1);
        return grid;
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
