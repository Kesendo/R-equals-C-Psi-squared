using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the running world (the diagonal protocol, ClaudeTasks/DIAGONAL_PROTOCOL_GAME.md).
// Step 1: the world splits. A field of possibilities loses, per tick, each pair's own Pair rate; the
// diagonal (k=0) stays (structure), the off-diagonal fades (novelty). Nothing interpreted, just the numbers.
public class FieldTests
{
    const double G = 0.5;
    static readonly World W = new();

    // the one question (rule 2): a pair fades by exactly its Pair rate (1 + Re lambda * dt); k=0 stays.
    [Fact]
    public void Field_Step_Fades_Each_Pair_By_Its_Pair_Rate()
    {
        var f = new Field(W, 2, G);
        int i = 0b01, j = 0b10;                         // popcount(01 ^ 10) = 2, a novelty pair
        f[i, j] = 1.0;
        f[0, 0] = 1.0;                                  // a structure pair, disagreement 0
        double dt = 0.05;

        f.Step(dt);

        double expected = 1.0 + new Pair(W, i, j, G).Rate * dt;   // Re lambda = -2 g k
        Assert.Equal(expected, f[i, j], 12);
        Assert.Equal(1.0, f[0, 0], 12);                 // structure costs nothing, untouched
    }

    // the split (rules 1-3): under the watching, structure holds and novelty is culled.
    [Fact]
    public void Field_Step_Keeps_Structure_Culls_Novelty()
    {
        var f = new Field(W, 2, G);
        f.SeedUniform();
        double structure0 = f.Structure, novelty0 = f.Novelty;

        f.Step(0.05);

        Assert.Equal(structure0, f.Structure, 12);      // the diagonal (k=0) never decays
        Assert.True(f.Novelty < novelty0);              // the off-diagonal is eroded
    }

    // the ontology: the field owns its split (left), and inherits the frame from the World (right).
    [Fact]
    public void Field_Owns_The_Split_And_Inherits_The_Frame()
    {
        var f = new Field(W, 2, G);
        Assert.Equal(new[] { "structure", "novelty" }, f.Own);
        Assert.Equal(new[] { "x", "y", "z" }, f.Inherited);
    }

    // the census by disagreement at the seed is the adopted bare sectors (Grading A palindrome [4,8,4]).
    [Fact]
    public void Field_WeightByDisagreement_Matches_Bare_Sectors_At_Uniform_Seed()
    {
        var f = new Field(W, 2, G);
        f.SeedUniform();
        var byK = f.WeightByDisagreement();
        var bare = Redistribution.Bare(2);              // [4, 8, 4]
        for (int k = 0; k <= 2; k++)
            Assert.Equal(bare[k], byK[k], 12);
    }

    // the world runs forward: each Step advances the field's own clock by dt.
    [Fact]
    public void Field_Step_Advances_The_Clock()
    {
        var f = new Field(W, 2, G);
        Assert.Equal(0.0, f.T, 12);
        f.Step(0.05);
        Assert.Equal(0.05, f.T, 12);
        f.Step(0.05);
        Assert.Equal(0.10, f.T, 12);
    }

    // the knower's cut + the mirror: the immortal diagonal (k=0) is held; of the off-diagonal, only the
    // upper triangle is stepped -- the lower is its transpose twin (rho = rho-dagger), derived for free.
    [Fact]
    public void Field_Steps_Only_The_Upper_Triangle_Mirroring_The_Rest()
    {
        var f = new Field(W, 2, G);
        var bare = Redistribution.Bare(2);              // [4, 8, 4]
        int offDiagonal = bare.Skip(1).Sum();           // 12 = 4^N - 2^N
        Assert.Equal(bare[0], f.ImmortalCount);         // 4 = 2^N, the held diagonal (k=0)
        Assert.Equal(offDiagonal / 2, f.AliveCount);    // 6, the upper triangle; the mirror gives the other 6
    }

    // the law that justifies the cut: reading either twin gives the same weight, every tick (rho = rho-dagger).
    [Fact]
    public void Field_Lower_Triangle_Mirrors_The_Upper()
    {
        var f = new Field(W, 2, G);
        f.SeedUniform();
        f.Step(0.05);
        f.Step(0.05);
        for (int i = 0; i < f.Dim; i++)
            for (int j = 0; j < f.Dim; j++)
                Assert.Equal(f[i, j], f[j, i], 12);
    }
}
