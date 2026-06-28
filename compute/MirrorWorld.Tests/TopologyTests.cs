using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the geometry (who can reach whom). The handshake rides these bonds; every one is
// excitation-conserving, so they share the block structure (F63) -- the geometry shapes how novelty
// propagates, not which cells are forbidden. Nothing interpreted, just the bond counts.
public class TopologyTests
{
    [Fact]
    public void Bond_Counts_Per_Geometry()
    {
        Assert.Equal(3, Topology.Chain(4).Length);       // N-1, a line
        Assert.Equal(4, Topology.Ring(4).Length);        // N, the line closed
        Assert.Equal(3, Topology.Star(4).Length);        // N-1, hub to every spoke
        Assert.Equal(6, Topology.Complete(4).Length);    // N(N-1)/2, all-to-all
    }

    [Fact]
    public void Ring_Is_The_Chain_For_Two()
    {
        Assert.Equal(Topology.Chain(2).Length, Topology.Ring(2).Length);   // a 2-ring is just the one bond
    }
}
