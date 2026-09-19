namespace MirrorWorld;

// A two-hand model for e^(lambda t) = e^(-alpha t) * e^(i omega t): the radial
// hand records decay alpha (scaled here by Gamma), and the angular hand records
// rotation omega (scaled here by J). Q=J/Gamma is a dimensionless dial comparing those hands
// when Gamma is nonzero; ThetaDeg=atan(Q) reports the dial angle.
// Gamma=0 is the no-decay limit of this clock and is outside F95's finite positive-b domain.
// The dial supplies no F95 ancestry and no quarter/half anchor.
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
