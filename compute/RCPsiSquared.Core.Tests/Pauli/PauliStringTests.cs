using System.Numerics;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Tests.Pauli;

public class PauliStringTests
{
    [Fact]
    public void TwoSiteStrings_HaveExpectedShape()
    {
        // I⊗X is 4×4 with X in upper-left 2×2 and X in lower-right.
        var ix = PauliString.Build(new[] { PauliLetter.I, PauliLetter.X });
        Assert.Equal(4, ix.RowCount);
        Assert.Equal(4, ix.ColumnCount);
        // (I⊗X)|01⟩ = |00⟩, (I⊗X)|00⟩ = |01⟩, etc. → entries at (0,1), (1,0), (2,3), (3,2)
        Assert.Equal(1.0, ix[0, 1].Real, 12);
        Assert.Equal(1.0, ix[1, 0].Real, 12);
        Assert.Equal(1.0, ix[2, 3].Real, 12);
        Assert.Equal(1.0, ix[3, 2].Real, 12);
    }

    [Fact]
    public void SiteOp_X_AtSite0_OfThreeQubits_HasCorrectShape()
    {
        var op = PauliString.SiteOp(N: 3, site: 0, PauliLetter.X);
        Assert.Equal(8, op.RowCount);
        // X⊗I⊗I in big-endian: site 0 is most significant bit, so flips bit 2 of the 3-bit index.
        // |000⟩ → |100⟩: index 0 → 4, etc.
        Assert.Equal(1.0, op[4, 0].Real, 12);
        Assert.Equal(1.0, op[0, 4].Real, 12);
    }

    [Fact]
    public void PauliStringTraces_ObeyOrthogonalityRelation()
    {
        // Tr(σ_α^† σ_β) = 2^N δ_αβ for two-letter Pauli strings.
        int N = 2;
        int d = 1 << N;
        long d2 = 1L << (2 * N);
        for (long a = 0; a < d2; a++)
        {
            for (long b = 0; b < d2; b++)
            {
                var sa = PauliString.Build(PauliIndex.FromFlat(a, N));
                var sb = PauliString.Build(PauliIndex.FromFlat(b, N));
                var trace = (sa * sb).Trace();
                double expected = (a == b) ? d : 0;
                Assert.Equal(expected, trace.Real, 10);
                Assert.Equal(0.0, trace.Imaginary, 10);
            }
        }
    }
}
