using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the geometry (who can reach whom). The handshake rides these bonds; every one is
// excitation-conserving, so they share the block structure (F63) -- the geometry shapes how novelty
// propagates, not which cells are forbidden. Check the actual undirected edges, including uniqueness.
public class TopologyTests
{
    static (int a, int b)[] Canonical((int a, int b)[] bonds) => bonds
        .Select(edge => (a: Math.Min(edge.a, edge.b), b: Math.Max(edge.a, edge.b)))
        .OrderBy(edge => edge.a).ThenBy(edge => edge.b).ToArray();

    [Fact]
    public void Bond_Sets_Per_Geometry()
    {
        Assert.Equal(new[] { (0, 1), (1, 2), (2, 3) }, Canonical(Topology.Chain(4)));
        Assert.Equal(new[] { (0, 1), (0, 3), (1, 2), (2, 3) }, Canonical(Topology.Ring(4)));
        Assert.Equal(new[] { (0, 1), (0, 2), (0, 3) }, Canonical(Topology.Star(4)));
        Assert.Equal(new[] { (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3) },
            Canonical(Topology.Complete(4)));
    }

    [Fact]
    public void Ring_Is_The_Chain_For_Two()
    {
        Assert.Equal(new[] { (0, 1) }, Canonical(Topology.Ring(2)));       // no doubled bond
    }

    [Fact]
    public void Named_RejectsUnknownGeometryInsteadOfComputingAChain()
    {
        Assert.Equal(Topology.Chain(4), Topology.Named("chain", 4));
        Assert.Throws<ArgumentException>(() => Topology.Named("strar", 4));
    }
}
