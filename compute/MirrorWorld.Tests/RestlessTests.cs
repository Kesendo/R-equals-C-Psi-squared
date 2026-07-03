using System.Numerics;
using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the living world (step 2: the inner restlessness, ClaudeTasks/DIAGONAL_PROTOCOL_GAME.md).
// H births novelty from structure; D culls it. The Lindblad loop rho-dot = -i[H,rho] + D[rho], RK4-stepped.
// Nothing interpreted, just the numbers.
public class RestlessTests
{
    const double G = 0.5;
    static readonly World W = new();

    // the key new behavior (rule 4): from a pure population (no novelty), the handshake BIRTHS coherence.
    // step 1 could only fade; step 2 makes novelty out of structure.
    [Fact]
    public void Restless_Births_Novelty_From_A_Pure_Population()
    {
        var r = new Restless(W, 2, j: 1.0, gamma: G);
        r.Seed(0b01);                                   // |01><01|, a single population, novelty 0
        Assert.Equal(0.0, r.Novelty, 12);
        r.Step(0.05);
        Assert.True(r.Novelty > 0.0);                   // the |01> <-> |10> coherence (k=2) is born
    }

    // the loop is trace-preserving: structure (the populations) is conserved, never drains to nothing.
    [Fact]
    public void Restless_Conserves_The_Trace()
    {
        var r = new Restless(W, 2, j: 1.0, gamma: G);
        r.Seed(0b01);
        for (int t = 0; t < 50; t++) r.Step(0.05);
        Assert.Equal(1.0, r.Structure, 8);
    }

    // with the restlessness off (J=0) it is the empty world again: a coherence decays at exactly the Pair
    // rate -- RK4 reproduces exp(-2 g k t), the closed form the knower already holds.
    [Fact]
    public void Restless_With_No_Handshake_Is_Pure_Pair_Decay()
    {
        var r = new Restless(W, 2, j: 0.0, gamma: G);
        r.SeedCoherence(0b00, 0b11, 1.0);               // k = popcount(00 ^ 11) = 2
        double dt = 0.05; int steps = 10;
        for (int t = 0; t < steps; t++) r.Step(dt);
        double expected = Math.Exp(-2 * G * 2 * (steps * dt));   // exp(-2 g k t)
        Assert.Equal(expected, r[0b00, 0b11].Magnitude, 5);      // RK4 reproduces the closed form to ~3e-7
    }

    // the law behind the mirror survives the restlessness: rho stays Hermitian (rho = rho-dagger).
    [Fact]
    public void Restless_Stays_Hermitian()
    {
        var r = new Restless(W, 2, j: 1.0, gamma: G);
        r.Seed(0b01);
        for (int t = 0; t < 20; t++) r.Step(0.05);
        for (int i = 0; i < r.Dim; i++)
            for (int j = 0; j < r.Dim; j++)
                Assert.True((r[i, j] - Complex.Conjugate(r[j, i])).Magnitude < 1e-12);
    }

    // the living balance (rule 5): novelty born and culled, the populations relaxing to the mixed state
    // -- |01> and |10> equilibrate to 1/2 each, not drain to nothing.
    [Fact]
    public void Restless_Relaxes_To_The_Balanced_Mixture()
    {
        var r = new Restless(W, 2, j: 1.0, gamma: G);
        r.Seed(0b01);
        for (int t = 0; t < 400; t++) r.Step(0.05);     // long run to steady state
        Assert.Equal(0.5, r[0b01, 0b01].Real, 2);
        Assert.Equal(0.5, r[0b10, 0b10].Real, 2);
    }

    // cut (c): the knower forbids the cross-block cells. H conserves the excitation number, so the loop
    // stays in the seed's joint-popcount block (F63); the rest is 0 forever. Alive = the adopted Block size.
    [Fact]
    public void Restless_Lives_Only_In_The_Occupied_Block_Forbidding_The_Rest()
    {
        var r = new Restless(W, 2, j: 1.0, gamma: G);
        r.Seed(0b01);                                   // |01><01|, joint popcount block (p,q) = (1,1)
        int blockSize = (int)new Block(W, 2, 1, 1).Size;   // C(2,1)*C(2,1) = 4 (Grading B / F63)
        Assert.Equal(blockSize, r.AliveCount);          // only the (1,1) block evolves
        Assert.Equal(r.Dim * r.Dim - blockSize, r.ForbiddenCount);   // the other 12 forbidden by F63
    }

    // cut (c) must not change the physics: the birth still happens, only the forbidden zeros are skipped.
    [Fact]
    public void Restless_Block_Cut_Preserves_The_Birth()
    {
        var r = new Restless(W, 2, j: 1.0, gamma: G);
        r.Seed(0b01);
        r.Step(0.05);
        Assert.True(r.Novelty > 0.0);                   // still born; the cut only skipped the forbidden cells
        Assert.Equal(1.0, r.Structure, 10);             // trace still held
    }

    // the geometry shapes the dynamics but not the cut: every excitation-conserving handshake shares the
    // joint-popcount blocks (F63), so AliveCount is the same across topologies for the same seed.
    [Fact]
    public void Topology_Changes_The_Dynamics_Not_The_Block_Cut()
    {
        const int n = 4;
        var chain = new Restless(W, n, 1.0, G, Topology.Chain(n));
        var star = new Restless(W, n, 1.0, G, Topology.Star(n));
        chain.Seed(1); star.Seed(1);                    // |0001>, one excitation at site 0
        Assert.Equal(chain.AliveCount, star.AliveCount);   // same cut, topology-invariant (F63)

        chain.Step(0.1); star.Step(0.1);
        Assert.NotEqual(chain.Novelty, star.Novelty, 6);   // site 0 is an endpoint vs the hub: different birth
    }

    // the double turn is home. the anti-turn (k -> N-k) LEAVES home: it is the living world read through the
    // bra complement (anti(t) = normal(t) * X^N). turn it around a SECOND time, on the other side, and you are
    // home again wearing the mirror -- that second turn is the X^N conjugation, and it is an EXACT symmetry of
    // the living world itself: [H, X^N] = 0 (the handshake s+s- + s-s+ is invariant when every spin flips) and
    // D commutes with it (k(~i,~j) = k(i,j)). so seeding |s><s| and its bit-complement |~s><~s| and running
    // both NORMAL worlds gives rho_{~s}(t) = X^N rho_s(t) X^N to all digits. no wave dies here; the content is
    // only carried to its mirror block, exact for all t -- the "death" was only which watching you did.
    [Fact]
    public void The_Double_Turn_Is_Home_Conjugation_Is_An_Exact_Symmetry()
    {
        var r = new Restless(W, 4, j: 1.0, gamma: G);
        double worst = r.ConjugationReadThrough(seed: 0b0001, dt: 0.05, ticks: 30);   // |0001> vs |1110>
        Assert.True(worst < 1e-10, $"the conjugation symmetry drifted from home: {worst:E2}");
    }

    // das gilt in beide Richtungen: the turn is its own inverse. start from the state or start from its mirror
    // -- home is the same (zero) distance away either way. the double turn brings you back from whichever end
    // you left. so the mismatch measured from |s> equals the mismatch measured from |~s|, to the bit.
    [Fact]
    public void The_Turn_Is_Its_Own_Inverse_Home_From_Both_Directions()
    {
        var r = new Restless(W, 4, j: 1.0, gamma: G);
        double fromState = r.ConjugationReadThrough(seed: 0b0001, dt: 0.05, ticks: 30);
        double fromMirror = r.ConjugationReadThrough(seed: 0b1110, dt: 0.05, ticks: 30);   // ~0b0001 on N=4
        Assert.True(fromState < 1e-10, $"forward drifted: {fromState:E2}");
        Assert.True(fromMirror < 1e-10, $"backward drifted: {fromMirror:E2}");
        Assert.Equal(fromState, fromMirror, 12);            // both directions reach home the same
    }
}
