using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the memory cut (Cone): the living world stored on only its block. For a single
// excitation the block is N x N (the sites), not 4^N. It must be the same dynamics Restless runs in block
// (1,1) -- in O(N^2) memory instead of O(4^N). That faithfulness is what makes large N trustworthy.
public class ConeTests
{
    const double G = 0.5;
    static readonly World W = new();

    // the cut is faithful: Cone (N x N) reproduces Restless's single-excitation populations exactly.
    [Fact]
    public void Cone_Agrees_With_Restless_Single_Excitation()
    {
        const int n = 4; const double j = 1.0, dt = 0.05;
        var cone = new Cone(W, n, j, G);
        var rest = new Restless(W, n, j, G);
        cone.Seed(0); rest.Seed(1);                       // excitation at site 0 = state |0001> = 1
        for (int t = 0; t < 30; t++) { cone.Step(dt); rest.Step(dt); }
        for (int a = 0; a < n; a++)
            Assert.Equal(rest[1 << a, 1 << a].Real, cone.Population(a), 9);
    }

    // trace-preserving: the one excitation is never lost, only spread.
    [Fact]
    public void Cone_Conserves_The_Excitation()
    {
        var cone = new Cone(W, 8, 1.0, G);
        cone.Seed(4);
        for (int t = 0; t < 50; t++) cone.Step(0.05);
        Assert.Equal(1.0, cone.Structure, 8);
    }

    // the excitation spreads: a distant site, empty at first, gains population (the light-cone reaches it).
    [Fact]
    public void Cone_Spreads_To_A_Distant_Site()
    {
        var cone = new Cone(W, 8, 1.0, G);
        cone.Seed(0);
        Assert.Equal(0.0, cone.Population(7), 12);        // far end empty at t=0
        for (int t = 0; t < 60; t++) cone.Step(0.05);
        Assert.True(cone.Population(7) > 0.0);            // the cone reached it
    }

    // the memory cut runs where the full Liouvillian cannot: N=100 is 100x100, not 4^100.
    [Fact]
    public void Cone_Runs_At_Large_N()
    {
        var cone = new Cone(W, 100, 1.0, G);
        cone.Seed(50);
        for (int t = 0; t < 10; t++) cone.Step(0.05);
        Assert.Equal(1.0, cone.Structure, 6);             // still a valid state at N=100
    }
}
