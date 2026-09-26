using System.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Move 4: the (SE,DE) diabolics pair across the (SE,DE)↔(SE,w_{N−2}) cross-block fold, because that
/// fold is an EXACT antiunitary similarity. The branch-locus palindrome's bra bit-flip ρ[a,b]→ρ[a,b̄] maps the
/// (w1,w2) block to the (w1,N−2) block; these tests pin the matrix identity
/// L(1,N−2)(q̄) = −P·conj(L(1,2)(q))·Pᵀ − 2N·I to machine zero (so the whole Jordan structure is preserved across
/// the fold), reproduce the N=7 real-q diabolic pairing at its certified location, and read an N=7 Δ = 0.02 EP2
/// and its partner turning defective in lockstep. The witness for the cross-fold section of
/// experiments/F89_PATH_K_DIABOLIC.md.</summary>
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
        // The N=7 real-q diabolic (q* = 1.126448513252489811, λ* = −4.941857490410003, certified at Δ=0 by
        // XxzDeltaFlipTests.N7_RealQCrossing_SplitsIntoTwoEp2s_SeedFollowing) maps to the (1,5) partner at the fold
        // image −λ*−2N. Error model at q*: each gap is 0.8632·|δq| (|δq| ≤ u·q* from rounding q*) plus the
        // eigensolver's floor at a semisimple crossing, about 1e-13 here; 1e-9 is what a q* wrong in its ninth digit
        // would give. The two gaps are EVDs of exactly similar matrices, so they differ by at most twice that floor.
        var w = new CrossFoldSimilarityWitness();
        var (g12, gp, partnerLam) = w.PairedGapsAcrossTheFold(7, CrossFoldSimilarityWitness.N7CrossingQ,
            CrossFoldSimilarityWitness.N7CrossingLambda);
        Assert.Equal(-9.058142509589997, partnerLam, 12);
        Assert.True(g12 < 1e-9, $"(1,2) is not a coalescence at q*: gap {g12:E2}");
        Assert.True(gp < 1e-9, $"partner (1,5) is not a coalescence at the fold image: gap {gp:E2}");
        Assert.True(System.Math.Abs(g12 - gp) < 1e-11, $"the paired gaps differ: {g12:E2} vs {gp:E2}");

        // Off q* both open linearly with one slope (a crossing, not a square-root branch), 0.863183 per unit q;
        // the next Taylor term is O(δq) relative, so 1e-3 bounds it at δq ≤ 1e-4.
        foreach (var dq in new[] { 1e-4, 1e-5, 1e-6 })
        {
            var (a, b, _) = w.PairedGapsAcrossTheFold(7, CrossFoldSimilarityWitness.N7CrossingQ + dq,
                CrossFoldSimilarityWitness.N7CrossingLambda);
            Assert.True(System.Math.Abs(a / dq - 0.863183) < 1e-3, $"slope {a / dq:R} at dq={dq}");
            Assert.True(System.Math.Abs(a - b) < 1e-11, $"paired gaps differ off q*: {a:E3} vs {b:E3}");
        }
    }

    [Fact]
    public void N7_DeltaEp2_AndItsFoldPartner_TurnDefectiveInLockstep()
    {
        // Δ = 0.02 splits the N=7 crossing into two Jordan EP2s; the first, read on the full (1,2) block, and its
        // partner on the (1,5) block at q̄ carry the same multiplicities and the same departure (the similarity
        // leaves all three unchanged). Departure model: the cancellation in ‖A‖² − Σ|λ|² costs u·‖A‖²/dep ≈ 4e-12
        // at dep 0.0015, so 1e-9 separates rounding from a different EP.
        var (ep, _) = XxzCoherenceBlock.CertifySplitUnderDelta(7, new Complex(1.1264, 0), new Complex(-4.942, 0), 0.02);
        Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Defective, ep.Verdict);
        var (src, par) = new CrossFoldSimilarityWitness().FoldedCharacters(7, ep.QCandidate, ep.LambdaCandidate, 0.02);
        Assert.Equal(2, src.Algebraic);
        Assert.Equal(1, src.Geometric);
        Assert.Equal(2, par.Algebraic);
        Assert.Equal(1, par.Geometric);
        Assert.True(System.Math.Abs(src.Departure - par.Departure) < 1e-9, $"{src.Departure:R} vs {par.Departure:R}");
        Assert.True(System.Math.Abs(src.Departure - ep.Departure) < 1e-9, $"full block {src.Departure:R} vs R-even sector {ep.Departure:R}");
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
        // integrability-independent. Where Δ splits a diabolic into EP2s, the partners stay similar, so they split alike.
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
