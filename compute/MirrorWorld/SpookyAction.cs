namespace MirrorWorld;

// Spooky action at a distance, composed from the atoms -- the k=2 sighting, the middle of the triptych
// ("k=1 is the pattern, k=2 the pair, k=N the cat", docs/quantum/SPOOKY_ACTION_TRANSLATED.md). It IS
// Field at N=2 holding the Bell-pair skeleton: the two definite records |00><00| + |11><11| on the
// diagonal (the sock drawer, immortal) and the one coherence |00><11| between them (the spook,
// disagreement k=2) paying -2*gamma*2 = -2*(gamma_A + gamma_B) at uniform gamma -- everything the pair
// holds beyond a classical pair of socks rides on that single entry. The one NEW read over Cat at N=2
// is the two Marginals, the local pages: a k=2 disagreement is invisible to any one-site page (F70
// kinematics -- and the Bell skeleton carries no k=1 content), so each page holds its immortal record
// while the spook between them decays. Nothing else is computed; the MEANING (no action -- the excess
// is what the between holds beyond the two pages) lives in SPOOKY_ACTION_TRANSLATED.md, not here.
public sealed class SpookyAction : GameObject
{
    public double Gamma { get; }
    readonly Field cloud;
    readonly Pair carrier;

    public SpookyAction(World world, double gamma) : base(world)
    {
        Gamma = gamma;
        cloud = new Field(world, 2, gamma);
        cloud[0b00, 0b00] = 1.0;                       // |00><00|: one definite record (bare weights, as DoubleSlit/Cat)
        cloud[0b11, 0b11] = 1.0;                       // |11><11|: the other record
        cloud[0b00, 0b11] = 1.0;                       // |00><11|: the spook, disagreement k=2
        carrier = new Pair(world, 0b00, 0b11, gamma);
        PageA = new Marginal(cloud, new[] { 0 });
        PageB = new Marginal(cloud, new[] { 1 });
    }

    public void Watch(double dt) => cloud.Step(dt);    // one tick of the watching
    public double T => cloud.T;

    // left: the phenomenon read off the composed atoms.
    public double SockDrawer => cloud.Structure;       // the classical 00/11 records: immortal
    public double Spook => cloud.Novelty;              // the k=2 coherence: everything beyond the socks
    public double SpookRate => carrier.Rate;           // -2*gamma*2 = -2*(gamma_A + gamma_B) at uniform gamma
    public Marginal PageA { get; }                     // one local page: blind to the spook (F70)
    public Marginal PageB { get; }                     // the other page: blind the same way

    public override IReadOnlyList<string> Own => new[] { "sockDrawer", "spook", "pages" };
}
