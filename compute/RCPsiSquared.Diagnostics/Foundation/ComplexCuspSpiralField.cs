using System.Globalization;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The Object Manager's telescope onto the interior axis in 2D: the complex-plane spiral into
/// the cusp circle |CΨ| = ¼. Where <see cref="InteriorHorizonField"/> reads the approach on the real
/// line, this reads it in the plane, when a common Z-drift Ω winds the coherence. The cusp point becomes
/// a circle; every spiral crosses it at the same time (the radial law is Ω-independent), and only the
/// crossing angle is free, the one IBM Kingston steered on demand.
///
/// <para>Five readings for one (γ, Ω, φ₀): the cusp circle (the contract in 2D, the point seen
/// edge-on); the spiral itself (CΨ_com(t) winding in, crossing the circle); the winding read along a
/// geometric Ω-ladder (the crossing time flat, the crossing angle moving, the steerable freedom); the
/// hardware (the two Kingston spirals and the on-demand steering, from the Confirmations registry); and
/// the F95 √-kinship (the angular winding carrying the same √-form the 1D heading does, kept at the
/// label). N-free, closed-form. The horizon is a dephasing fold, never gravitational.</para></summary>
public sealed class ComplexCuspSpiralField : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private readonly double _gamma;
    private readonly double _omega;
    private readonly double _phi0;
    private readonly double[] _omegaLadder;   // geometric, ascending
    private readonly double _tMaxFactor;      // the trajectory runs to tMaxFactor × the crossing time

    public ComplexCuspSpiralField(double gamma = 0.5, double omega = 0.4, double phi0 = 0.0,
        int omegaPoints = 9, double tMaxFactor = 4.0)
    {
        if (gamma <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), $"gamma must be positive; got {gamma}");
        if (omega < 0) throw new ArgumentOutOfRangeException(nameof(omega), $"omega must be non-negative; got {omega}");
        if (omegaPoints < 2) throw new ArgumentOutOfRangeException(nameof(omegaPoints), $"need at least two Ω points; got {omegaPoints}");
        if (tMaxFactor <= 1.0) throw new ArgumentOutOfRangeException(nameof(tMaxFactor), $"tMaxFactor must exceed 1 (run past the crossing); got {tMaxFactor}");
        _gamma = gamma;
        _omega = omega;
        _phi0 = phi0;
        _tMaxFactor = tMaxFactor;
        // The ladder spans a decade of winding up to the chosen Ω (a fixed band when Ω = 0), so the
        // crossing-angle-vs-Ω reading has range. Geometric in Ω.
        double hi = omega > 0 ? omega : 1.0;
        _omegaLadder = GeometricLadder(hi / 10.0, hi, omegaPoints);
    }

    private static double[] GeometricLadder(double lo, double hi, int count)
    {
        var grid = new double[count];
        double logLo = Math.Log(lo), logHi = Math.Log(hi);
        for (int i = 0; i < count; i++)
            grid[i] = Math.Exp(logLo + (logHi - logLo) * i / (count - 1));
        return grid;
    }

    private double[] TimeGrid()
    {
        double tMax = _tMaxFactor * ComplexCuspSpiral.CrossingTime(_gamma);
        const int points = 200;
        var grid = new double[points];
        for (int i = 0; i < points; i++) grid[i] = tMax * i / (points - 1);
        return grid;
    }

    public string DisplayName =>
        $"ComplexCuspSpiralField (the interior axis in 2D, γ={_gamma.ToString("0.###", Inv)}, Ω={_omega.ToString("0.###", Inv)}, winding Ω/4γ={ComplexCuspSpiral.WindingRate(_gamma, _omega).ToString("0.##", Inv)})";

    public string Summary
    {
        get
        {
            double tc = ComplexCuspSpiral.CrossingTime(_gamma);
            double angle = ComplexCuspSpiral.CrossingArgument(_gamma, _omega, _phi0) * 180.0 / Math.PI;
            return $"the cusp ¼ is a circle |CΨ|=¼; the spiral crosses it at t={tc.ToString("0.###", Inv)} (Ω-independent) " +
                   $"at angle {angle.ToString("0.#", Inv)}° (Ω-set). Every spiral crosses the same circle; only the angle is free, " +
                   "the one Kingston steered. The 1D point seen in 2D; no gravity.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // 1. The cusp circle (the contract in 2D): |CΨ| = ¼, the point seen edge-on.
            const int circlePts = 200;
            var cx = new double[circlePts];
            var cy = new double[circlePts];
            for (int i = 0; i < circlePts; i++)
            {
                double a = 2.0 * Math.PI * i / (circlePts - 1);
                cx[i] = ComplexCuspSpiral.CircleRadius * Math.Cos(a);
                cy[i] = ComplexCuspSpiral.CircleRadius * Math.Sin(a);
            }
            yield return new InspectableNode(
                displayName: "the cusp circle (the contract in 2D)",
                summary: $"|CΨ|=¼ is a circle (radius {ComplexCuspSpiral.CircleRadius}, center 0): the 1D cusp point seen edge-on. Inert, universal; every spiral crosses it.",
                payload: new InspectablePayload.Curve("the ¼-circle", cx, cy, "Re(CΨ)", "Im(CΨ)"));

            // 2. The spiral (one trajectory): CΨ_com(t) winding into the circle.
            var tGrid = TimeGrid();
            var re = tGrid.Select(t => ComplexCuspSpiral.Re(_gamma, _omega, _phi0, t)).ToArray();
            var im = tGrid.Select(t => ComplexCuspSpiral.Im(_gamma, _omega, _phi0, t)).ToArray();
            double turns = ComplexCuspSpiral.WindingNumber(_omega, tGrid[^1]);
            yield return new InspectableNode(
                displayName: "the spiral (one trajectory, winding in)",
                summary: $"CΨ_com(t) = |CΨ|·e^(i(φ₀−Ωt)) from 1/3 inward, {turns.ToString("0.##", Inv)} turns over the grid; crosses the ¼-circle once. Ω=0 would be the real-axis line (the 1D axis).",
                payload: new InspectablePayload.Curve("the spiral CΨ_com(t)", re, im, "Re(CΨ)", "Im(CΨ)"));

            // 3. The winding (geometric Ω-ladder): the crossing time flat, the crossing angle moving.
            double tCrossFlat = ComplexCuspSpiral.CrossingTime(_gamma);
            var angles = _omegaLadder.Select(w => ComplexCuspSpiral.CrossingArgument(_gamma, w, _phi0) * 180.0 / Math.PI).ToArray();
            yield return new InspectableNode(
                displayName: "the winding (the angle is the free thing)",
                summary: $"across Ω the crossing time is flat at t={tCrossFlat.ToString("0.###", Inv)} (the radial law is Ω-independent); the crossing angle sweeps {angles[0].ToString("0.#", Inv)}°..{angles[^1].ToString("0.#", Inv)}°. The steerable freedom (Kingston f95_angle_steering).",
                payload: new InspectablePayload.Curve("crossing angle vs Ω", _omegaLadder, angles, "Ω", "crossing angle°"));

            // 4. The hardware (Kingston): the two observed spirals + the on-demand steering.
            yield return new InspectableNode(
                displayName: "the hardware (the Kingston spirals)",
                summary: "IBM Kingston 2026-04: Pair A spirals clockwise (arg −8°→−60°), Pair B counter-clockwise (+15°→+79°), both crossing |CΨ|=¼ " +
                         "(f25_cusp_trajectory, f57_kdwell_gamma_invariance). 2026-05: the crossing angle steered on demand by an injected Ω " +
                         "(f95_angle_steering_kingston_may2026, three crossings, 6.8°–15.7° residual). The angle is real and controllable.");

            // 5. The F95 √-kinship (the reading): the angular winding carries the 1D heading's √-form.
            yield return new InspectableNode(
                displayName: "the F95 √-kinship (the reading, at the label)",
                summary: "the radial dwell is F57, the angular winding carries the F95 √-form, and ¼ is the discriminant zero where they meet (CPSI_COMPLEX_PLANE.md). " +
                         "Sibling of the interior heading θ=arctan(√(4CΨ−1)); the solid cusp/EP F95 algebra is typed in TransitionBridgeF95SiblingClaim. Kept at the label, not asserted here.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
