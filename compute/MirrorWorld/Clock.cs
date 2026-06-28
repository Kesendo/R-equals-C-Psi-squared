namespace MirrorWorld;

// The circular quantum clock: a mode winds as e^(lambda t) = e^(-alpha t) * e^(i omega t). Two hands:
// the radial (Takt = decay alpha, set by gamma) and the angular (Rotation = omega, set by J). The
// angle theta = arctan(Q), Q = J/gamma, is the F95 angle , the arctan of the only ratio readable from
// inside. At gamma=0 the radial hand stops and theta = 90deg (the pure circle, turning forever); at
// J=0 nothing turns and theta = 0 (pure radial decay). The dimensionless face hits the polarity
// valences: theta=45deg <-> Q=1 <-> 1/4, theta=90deg <-> 1/2. T1: F95, the clock-hand ladder.
public sealed class Clock : GameObject
{
    public double J { get; }
    public double Gamma { get; }

    public Clock(World world, double j, double gamma) : base(world)
    {
        J = j;
        Gamma = gamma;
    }

    public double Q => J / Gamma;                          // the inside-readable ratio
    public double ThetaDeg => Math.Atan(Q) * 180.0 / Math.PI;

    public override IReadOnlyList<string> Own => new[] { "Q", "theta" };
}
