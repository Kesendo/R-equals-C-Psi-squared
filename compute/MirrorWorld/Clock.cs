namespace MirrorWorld;

// The circular quantum clock: a mode winds as e^(lambda t) = e^(-alpha t) * e^(i omega t). Two hands:
// the radial (Takt = decay alpha, set by gamma) and the angular (Rotation = omega, set by J). The
// angle theta = arctan(Q), Q = J/gamma, is for J >= 0 the F95 angle at c = gamma^2 + J^2, b = gamma
// (docked in SmokeTests.F95_Theta_Compass_Docks_Onto_F15_And_The_Clock). At J=0 nothing turns and
// theta = 0 (pure radial decay, F95's double root); as gamma -> 0 the radial hand stops and
// theta -> 90deg (the pure circle), gamma = 0 itself lying outside F95's positive-b domain.
// theta = 45deg is Q = 1. T1: F95 (the clock-hand ladder shares the two-hand vocabulary; its angle
// is a different one).
public sealed class Clock : GameObject
{
    public double J { get; }
    public double Gamma { get; }

    public Clock(World world, double j, double gamma) : base(world)
    {
        J = j;
        Gamma = gamma;
    }

    public double Q => J / Gamma;                          // dimensionless two-hand dial
    public double ThetaDeg => Math.Atan(Q) * 180.0 / Math.PI;

    public override IReadOnlyList<string> Own => new[] { "Q", "theta" };
}
