namespace MirrorWorld;

// The circular quantum clock: a mode winds as e^(lambda t) = e^(-alpha t) * e^(i omega t). Two hands:
// the radial (Takt = decay alpha, set by gamma) and the angular (Rotation = omega, set by J). Nobody here
// knows what time is, and the clock does not say; it holds the two hands and the angle between them.
//
// One mode (no N): the mode that turns at J while it fades at gamma, both its own rates. The angle
// theta = arctan(Q), Q = J/gamma, is for J >= 0 the F95 angle at c = gamma^2 + J^2, b = gamma (docked in
// SmokeTests.F95_Theta_Compass_Docks_Onto_F15_And_The_Clock). At J=0 nothing turns and theta = 0 (pure
// radial decay, F95's double root); as gamma -> 0 the radial hand stops and theta -> 90deg (the pure
// circle), gamma = 0 itself lying outside F95's positive-b domain. theta = 45deg is Q = 1.
//
// A chain (N given: the open XY chain Restless and Cone run, hop J, uniform Z light at the channel rate
// gamma, a Pair's rate being -2 gamma k): ClockHandLadderClaim adopted. The object owns the Takt hand
// 2 gamma (a Pair's rate at one disagreement), the coherence hand Formulas.OmegaMem (the band edge
// 2J cos(pi/(N+1)) for N >= 3, the pulled 2 sqrt(J^2 - gamma^2) at N = 2) and the angle between them,
// arctan(Q cos(pi/(N+1))) for N >= 3 and arctan(sqrt(Q^2 - 1)) at N = 2. Where these are the hands the
// chain shows is the claim's to say (ClockHandLadderClaim; the two-site floor in
// PROOF_CHAIN_GAP_DOMINANCE 4.2). HandsHeld reads the claim's sufficient side, Q > Formulas.Qstar(N);
// below it ThetaDeg is still the angle of the hands owned here. T1: F95, ClockHandLadderClaim.
public sealed class Clock : GameObject
{
    public double J { get; }
    public double Gamma { get; }
    public int? N { get; }

    public Clock(World world, double j, double gamma, int? n = null) : base(world)
    {
        if (n is < 2) throw new ArgumentOutOfRangeException(nameof(n), n, "a chain clock needs N >= 2");
        J = j;
        Gamma = gamma;
        N = n;
    }

    public double Q => J / Gamma;                          // dimensionless two-hand dial

    public double? TaktHand => N is null ? null : 2.0 * Gamma;                        // one unit of light
    public double? CoherenceHand => N is int n ? Formulas.OmegaMem(n, J, Gamma) : null;
    public bool? HandsHeld => N is int n ? Q > Formulas.Qstar(n) : null;

    // One mode: arctan(Q). A chain: the angle between its two hands.
    public double ThetaDeg => (N is null
        ? Math.Atan(Q)
        : Math.Atan2(CoherenceHand!.Value, TaktHand!.Value)) * 180.0 / Math.PI;

    public override IReadOnlyList<string> Own => N is null
        ? new[] { "Q", "theta" }
        : new[] { "Q", "theta", "takt", "coherence" };
}
