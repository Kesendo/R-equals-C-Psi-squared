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

    [Theory]
    [InlineData(5, 0.5)]
    [InlineData(6, 0.7)]
    [InlineData(7, 1.0)]
    [InlineData(9, -0.5)]
    public void CrossFold_SurvivesXxzAnisotropy_RealQ(int nBlock, double delta)
    {
        // The (q,Δ) extension: the antiunitary similarity holds for the FULL interacting XXZ block at Δ≠0, NOT
        // just the integrable XY one (the Δ·ZZ term is even under the global bit-flip). The fold is therefore
        // integrability-independent: it survives the very anisotropy that makes the diabolics defective.
        var r = new CrossFoldSimilarityWitness().Read(nBlock, new Complex(1.0, 0), delta);
        Assert.Equal(delta, r.Delta);
        Assert.True(r.SimilarityResidual < 1e-9,
            $"cross-fold broke under XXZ at N={nBlock}, Δ={delta}: residual {r.SimilarityResidual:E2}");
    }

    [Theory]
    [InlineData(6, 0.7)]
    [InlineData(8, 0.4)]
    public void CrossFold_SurvivesXxzAnisotropy_ComplexQ(int nBlock, double delta)
    {
        var r = new CrossFoldSimilarityWitness().Read(nBlock, new Complex(0.6407, 0.180), delta);
        Assert.True(r.SimilarityResidual < 1e-9,
            $"cross-fold broke under XXZ at complex q, N={nBlock}, Δ={delta}: residual {r.SimilarityResidual:E2}");
    }

    [Fact]
    public void DeltaZero_Reading_MatchesTheLegacyReadOverload()
    {
        // The Δ=0 reading of the new overload must equal the original Read(n,q) (the delegation contract).
        var w = new CrossFoldSimilarityWitness();
        var q = new Complex(0.9, -0.2);
        Assert.Equal(w.Read(7, q).SimilarityResidual, w.Read(7, q, 0.0).SimilarityResidual, 12);
    }

    [Fact]
    public void LongitudinalZField_BreaksTheFold_TheBitFlipParityDiscriminant()
    {
        // The complementary control: a bit-flip-ODD perturbation (a longitudinal Z-field) breaks the fold, so the
        // survival result above is not vacuous. The discriminant is bit-flip parity: even (ZZ) survives, odd
        // (field) breaks. Residual is O(1), not machine zero.
        double[] field = { 0.4, -0.3, 0.6, 0.2, -0.5, 0.1 };
        double res = new CrossFoldSimilarityWitness().ReadFieldControlResidual(6, new Complex(1.3, 0), field);
        Assert.True(res > 1.0, $"a longitudinal Z-field should break the cross-fold, but residual was only {res:E2}");
    }

    [Theory]
    [InlineData(6, 2, 3, 0.6)]    // bra leg at wKet=2 (F89d generalized past wKet=1)
    [InlineData(7, 2, 2, 0.0)]
    [InlineData(7, 3, 2, 1.0)]
    public void BraLeg_IsExact_AtGeneralKetWeight(int n, int wKet, int wBra, double delta)
    {
        double res = new CrossFoldSimilarityWitness().BraLegResidual(n, wKet, wBra, new Complex(1.3, -0.2), delta);
        Assert.True(res < 1e-9, $"bra-leg broke at N={n}, ({wKet},{wBra}), Δ={delta}: residual {res:E2}");
    }

    [Theory]
    [InlineData(6, 2, 3, 0.6)]    // the NEW ket leg (mirror of F89d on the ket index)
    [InlineData(7, 2, 2, 0.0)]
    [InlineData(7, 1, 3, 1.0)]
    public void KetLeg_IsExactAntiunitarySimilarity(int n, int wKet, int wBra, double delta)
    {
        double res = new CrossFoldSimilarityWitness().KetLegResidual(n, wKet, wBra, new Complex(1.3, -0.2), delta);
        Assert.True(res < 1e-9, $"ket-leg broke at N={n}, ({wKet},{wBra}), Δ={delta}: residual {res:E2}");
    }

    [Theory]
    [InlineData(6, 2, 3, 0.6)]    // the unitary global spin-flip QP = X^⊗N = Π²
    [InlineData(7, 2, 2, 0.5)]
    public void FullFlip_IsUnitarySpinFlipSimilarity(int n, int wKet, int wBra, double delta)
    {
        double res = new CrossFoldSimilarityWitness().FullFlipResidual(n, wKet, wBra, new Complex(1.3, -0.2), delta);
        Assert.True(res < 1e-9, $"full-flip (spin-flip) broke at N={n}, ({wKet},{wBra}), Δ={delta}: residual {res:E2}");
    }

    [Fact]
    public void Summary_StatesMove4Answered()
    {
        var s = new CrossFoldSimilarityWitness().Summary;
        Assert.Contains("similarity", s, System.StringComparison.OrdinalIgnoreCase);
        Assert.Contains("pair", s, System.StringComparison.OrdinalIgnoreCase);
        Assert.Contains("Δ", s);                                              // the Δ-robustness is stated
        Assert.Contains("leg", s, System.StringComparison.OrdinalIgnoreCase); // the two-leg Klein structure is stated
        Assert.False(string.IsNullOrWhiteSpace(s));
    }
}
