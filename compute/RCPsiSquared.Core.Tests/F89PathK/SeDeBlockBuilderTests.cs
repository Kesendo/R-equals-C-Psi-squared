using RCPsiSquared.Core.F89PathK;
using Xunit;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The shared (SE,DE) block primitive: the ONE builder of the (marked point, 2-subset)
/// Liouvillian sub-block (J=γ=1, diagonal −2 overlap / −6 no-overlap, hopping ±2 governed by the
/// graph), consumed by both TopologyGaloisWritabilityWitness and the Galois-vs-spectral-chaos witness.
/// Extracted so the block is built once, not re-implemented per consumer.</summary>
public class SeDeBlockBuilderTests
{
    [Theory]
    [InlineData(4, 24)]   // 4 · C(4,2) = 4·6
    [InlineData(5, 50)]   // 5 · C(5,2) = 5·10
    [InlineData(7, 147)]  // 7 · C(7,2) = 7·21
    public void Basis_HasDimension_NTimesNChoose2(int n, int expected)
    {
        Assert.Equal(expected, SeDeBlockBuilder.Basis(n).Count);
    }

    [Fact]
    public void Build_Diagonal_IsMinus2OnOverlapAndMinus6Otherwise()
    {
        // For each 2-subset {j,k} exactly 2 of the n marked points overlap (i=j or i=k).
        // n=5: C(5,2)=10 subsets ⟹ 20 overlap rows (−2), the remaining 30 no-overlap (−6).
        var (re, _) = SeDeBlockBuilder.Build(5, "complete");
        int minus2 = 0, minus6 = 0;
        int d = SeDeBlockBuilder.Basis(5).Count;
        for (int t = 0; t < d; t++)
        {
            if (re[t, t] == -2) minus2++;
            else if (re[t, t] == -6) minus6++;
        }
        Assert.Equal(20, minus2);
        Assert.Equal(30, minus6);
    }

    [Fact]
    public void Build_HoppingMatrix_IsSymmetric()
    {
        // im is the XY hopping (SE ket + DE bra), symmetric ⟹ L = re + i·im is complex symmetric.
        var (_, im) = SeDeBlockBuilder.Build(6, "chain");
        int d = im.GetLength(0);
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
                Assert.Equal(im[a, b], im[b, a]);
    }

    [Fact]
    public void Adjacency_Ring_IsACycle_EachSiteHasTwoNeighbours()
    {
        // the F89 ring is C_N: site 0 ~ {1, N−1}, not {1,2} (chain). This is the new topology the
        // shared builder adds over the witness's chain/star/complete.
        var adj = SeDeBlockBuilder.Adjacency(4, "ring");
        Assert.True(adj[0, 1]);
        Assert.True(adj[0, 3]);
        Assert.False(adj[0, 2]);
        for (int i = 0; i < 4; i++)
        {
            int deg = 0;
            for (int j = 0; j < 4; j++) if (adj[i, j]) deg++;
            Assert.Equal(2, deg);
        }
    }
}
