using System.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Move 4: the (SE,DE) diabolics pair across the (SE,DE)↔(SE,w_{N−2}) cross-block fold, because that
/// fold is an EXACT antiunitary similarity. The branch-locus palindrome's bra bit-flip ρ[a,b]→ρ[a,b̄] maps the
/// (w1,w2) block to the (w1,N−2) block; these tests pin the matrix identity
/// L(1,N−2)(q̄) = −P·conj(L(1,2)(q))·Pᵀ − 2N·I to machine zero (so the whole Jordan structure, hence every
/// diabolic's character, is preserved across the fold) and reproduce the N=7 real-q diabolic pairing. The
/// witness for the cross-fold section of experiments/F89_PATH_K_DIABOLIC.md.</summary>
public class CrossFoldSimilarityWitnessTests
{
    [Theory]
    [InlineData(4)]     // N=4: partner w_{N−2}=2=DE is the (SE,DE) block itself (the within-block self-fold)
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    public void CrossFold_IsExactAntiunitarySimilarity_RealQ(int nBlock)
    {
        var r = new CrossFoldSimilarityWitness().Read(nBlock, new Complex(1.0, 0));
        Assert.Equal(nBlock - 2, r.PartnerWBra);
        Assert.True(r.SimilarityResidual < 1e-9,
            $"cross-fold similarity not exact at N={nBlock}: residual {r.SimilarityResidual:E2}");
    }

    [Theory]
    [InlineData(5)]     // the identity is the F1 form L(1,N−2)(q̄) = −P conj(L(1,2)(q)) Pᵀ − 2N·I, so it holds
    [InlineData(7)]     // at COMPLEX q too (partner evaluated at the conjugate coupling)
    [InlineData(8)]
    public void CrossFold_IsExactAntiunitarySimilarity_ComplexQ(int nBlock)
    {
        var r = new CrossFoldSimilarityWitness().Read(nBlock, new Complex(0.6407, 0.180));
        Assert.True(r.SimilarityResidual < 1e-9,
            $"cross-fold similarity not exact at complex q, N={nBlock}: residual {r.SimilarityResidual:E2}");
    }

    [Fact]
    public void N7_RealQDiabolic_PairsAcrossTheFold()
    {
        // The N=7 real-q diabolic (q=1.1264, λ=−4.942) in (SE,DE)=(1,2) maps to the partner (1,5) at the fold
        // image −λ−2N = −9.058; the two coalescence gaps are equal (the exact similarity), and both small.
        var (g12, gp, partnerLam) = new CrossFoldSimilarityWitness().ReproducePairedDiabolic(7, 1.1264, -4.942);
        Assert.Equal(-9.058, partnerLam, 3);
        Assert.True(g12 < 1e-3, $"(1,2) is not a coalescence near λ=−4.942: gap {g12:E2}");
        Assert.True(gp < 1e-3, $"partner (1,5) is not a coalescence near −9.058: gap {gp:E2}");
        Assert.True(System.Math.Abs(g12 - gp) < 1e-9, $"the paired gaps differ: {g12:E2} vs {gp:E2}");
    }

    [Fact]
    public void Summary_StatesMove4Answered()
    {
        var s = new CrossFoldSimilarityWitness().Summary;
        Assert.Contains("similarity", s, System.StringComparison.OrdinalIgnoreCase);
        Assert.Contains("pair", s, System.StringComparison.OrdinalIgnoreCase);
        Assert.False(string.IsNullOrWhiteSpace(s));
    }
}
