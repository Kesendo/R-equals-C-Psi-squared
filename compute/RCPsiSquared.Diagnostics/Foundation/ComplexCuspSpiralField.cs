using System.Globalization;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The Object Manager's view of a named Bell+/pure-Z complex-coherence trajectory.
/// A common Z-drift Ω winds the phase while the F25 magnitude decays. The ring
/// |CΨ_com|=¼ is a radial readout set, not the period-one cardioid and not a recurrence root locus.
///
/// <para>Five readings for one (γ, Ω, φ₀): the radial ring; the spiral; an Ω-ladder (the crossing
/// time flat, the crossing argument moving, the steerable freedom); the Kingston records; and the
/// comparison with the recurrence and its F95 heading (experiments/CPSI_COMPLEX_PLANE.md). The closed
/// forms are N-free within this setup; the crossing is structural, never gravitational.</para></summary>
public sealed class ComplexCuspSpiralField : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private readonly double _gamma;
    private readonly double _omega;
    private readonly double _phi0;
    private readonly double[] _omegaLadder;   // geometric in |Ω|, ascending in magnitude, carrying Ω's sign
    private readonly double _tMaxFactor;      // the trajectory runs to tMaxFactor × the crossing time

    public ComplexCuspSpiralField(double gamma = 0.5, double omega = 0.4, double phi0 = 0.0,
        int omegaPoints = 9, double tMaxFactor = 4.0)
    {
        if (gamma <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), $"gamma must be positive; got {gamma}");
        if (omegaPoints < 2) throw new ArgumentOutOfRangeException(nameof(omegaPoints), $"need at least two Ω points; got {omegaPoints}");
        if (tMaxFactor <= 1.0) throw new ArgumentOutOfRangeException(nameof(tMaxFactor), $"tMaxFactor must exceed 1 (run past the crossing); got {tMaxFactor}");
        _gamma = gamma;
        _omega = omega;
        _phi0 = phi0;
        _tMaxFactor = tMaxFactor;
        // The ladder spans a decade of winding up to the chosen Ω (a fixed band when Ω = 0), so the
        // crossing-angle-vs-Ω reading has range. Geometric in |Ω|, carrying Ω's sign: the two
        // Kingston spirals this field advertises turn opposite ways (Pair A clockwise, Pair B
        // counter-clockwise), and with arg = φ₀ − Ω·t the counter-clockwise one needs Ω < 0.
        double sign = omega < 0 ? -1.0 : 1.0;
        double hi = omega != 0.0 ? Math.Abs(omega) : 1.0;
        var ladder = GeometricLadder(hi / 10.0, hi, omegaPoints);
        for (int i = 0; i < ladder.Length; i++) ladder[i] *= sign;
        _omegaLadder = ladder;
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
        $"ComplexCuspSpiralField (Bell+/pure-Z complex readout, γ={_gamma.ToString("0.###", Inv)}, Ω={_omega.ToString("0.###", Inv)}, winding Ω/4γ={ComplexCuspSpiral.WindingRate(_gamma, _omega).ToString("0.##", Inv)})";

    public string Summary
    {
        get
        {
            double tc = ComplexCuspSpiral.CrossingTime(_gamma);
            double angle = ComplexCuspSpiral.CrossingArgument(_gamma, _omega, _phi0) * 180.0 / Math.PI;
            return $"the Bell+/pure-Z magnitude reaches the radial ring |CΨ_com|=¼ at t={tc.ToString("0.###", Inv)} (Ω-independent) " +
                   $"with argument {angle.ToString("0.#", Inv)}° (Ω-set). The ring is a finite readout set; the recurrence boundary and cardioid are separate objects.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // 1. The selected radial readout ring |CΨ_com| = ¼.
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
                displayName: "the radial ring (finite readout set)",
                summary: $"|CΨ_com|=¼ is a radial ring (radius {ComplexCuspSpiral.CircleRadius}, center 0). It is not the F97 cardioid; only c=+¼ is that cardioid's cusp.",
                payload: new InspectablePayload.Curve("the quarter-radius ring", cx, cy, "Re(CΨ_com)", "Im(CΨ_com)"));

            // 2. The spiral (one trajectory): CΨ_com(t) winding into the circle.
            var tGrid = TimeGrid();
            var re = tGrid.Select(t => ComplexCuspSpiral.Re(_gamma, _omega, _phi0, t)).ToArray();
            var im = tGrid.Select(t => ComplexCuspSpiral.Im(_gamma, _omega, _phi0, t)).ToArray();
            double turns = ComplexCuspSpiral.WindingNumber(_omega, tGrid[^1]);
            yield return new InspectableNode(
                displayName: "the spiral (one trajectory, winding in)",
                summary: $"CΨ_com(t) = |CΨ_com|·e^(i(φ₀−Ωt)) from 1/3 inward, {turns.ToString("0.##", Inv)} turns over the grid; reaches the quarter-radius ring once. Ω=0 gives the ray at φ₀ (the real axis at φ₀=0, ideal Bell+).",
                payload: new InspectablePayload.Curve("the spiral CΨ_com(t)", re, im, "Re(CΨ)", "Im(CΨ)"));

            // 3. The winding (geometric Ω-ladder): the crossing time flat, the crossing angle moving.
            double tCrossFlat = ComplexCuspSpiral.CrossingTime(_gamma);
            var angles = _omegaLadder.Select(w => ComplexCuspSpiral.CrossingArgument(_gamma, w, _phi0) * 180.0 / Math.PI).ToArray();
            yield return new InspectableNode(
                displayName: "the winding (the angle is the free thing)",
                summary: $"across Ω the radial crossing time is flat at t={tCrossFlat.ToString("0.###", Inv)}; the argument sweeps {angles[0].ToString("0.#", Inv)}°..{angles[^1].ToString("0.#", Inv)}°. The steerable freedom (f95_angle_steering_kingston_may2026).",
                payload: new InspectablePayload.Curve("crossing angle vs Ω", _omegaLadder, angles, "Ω", "crossing angle°"));

            // 4. The hardware (Kingston): the two observed spirals + the on-demand steering.
            yield return new InspectableNode(
                displayName: "the hardware (the Kingston spirals)",
                summary: "IBM Kingston 2026-04: Pair A spirals clockwise (arg −8°→−60°), Pair B counter-clockwise (+15°→+79°), both reaching |CΨ_com|=¼ " +
                         "(f25_cusp_trajectory, f57_kdwell_gamma_invariance). 2026-05: the crossing angle steered on demand by an injected Ω " +
                         "(f95_angle_steering_kingston_may2026, three crossings, residuals 6.8° to 15.7° from one Lindblad+RZ model with in-situ γ). The argument is real and controllable.");

            // 5. The comparison: the recurrence and its F95 heading against this ring and its argument.
            yield return new InspectableNode(
                displayName: "the comparison (recurrence heading and steered argument)",
                summary: "the recurrence boundary c=¼ is the double root of z²−z+c=0, and its heading θ=arctan(√(4c−1)) is F95's θ(c; ½) on that quadratic. " +
                         "This radial ring is a different object, the Bell+/pure-Z magnitude reaching ¼, not a root locus. " +
                         "The steered argument arg(CΨ_com) shares F95's square-root form in its dwell arc length but is not θ itself (CPSI_COMPLEX_PLANE.md); " +
                         "the cusp/EP F95 algebra is typed in TransitionBridgeF95SiblingClaim.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
