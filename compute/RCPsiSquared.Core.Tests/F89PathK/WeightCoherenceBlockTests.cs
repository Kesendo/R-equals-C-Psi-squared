using System;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.SymmetryFamily;
using Xunit;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The general (wKet, wBra) coherence-block builder, its ZZ-bond sum, and the F89d cross-fold antiunitary
/// similarity EXTENDED to the interacting XXZ chain (Δ ≠ 0). The cross-fold identity
/// L(1,N−2)(q̄,Δ) = −P · conj(L(1,2)(q,Δ)) · Pᵀ − 2N·I holds to machine zero at every Δ, because the Δ·ZZ term is
/// EVEN under the global bit-flip (zz(b̄) = zz(b)), so the diabolic pairing is integrability-INDEPENDENT (it
/// survives the very anisotropy that kills the diabolics themselves). The discriminant is bit-flip parity: a
/// bit-flip-ODD perturbation (a longitudinal Z-field, fe(b̄) = −fe(b)) breaks the fold. The Core gate for
/// the cross-fold (q,Δ) extension; the live evidence is CrossFoldSimilarityWitness (Diagnostics).</summary>
public class WeightCoherenceBlockTests
{
    private static readonly int[] SweepN = { 4, 5, 6, 7, 8, 9 };

    // The cross-fold antiunitary-similarity residual at (q, Δ):
    // max over (t,u) of |L(1,N−2)(q̄,Δ)[Pt,Pu] − (−conj(L(1,2)(q,Δ)[t,u]) − 2N·δ)|. Zero ⟹ exact similarity.
    private static double CrossFoldResidual(int n, Complex q, double delta)
    {
        var l12 = WeightCoherenceBlock.Build(n, 1, 2, q, delta);
        var lpartner = WeightCoherenceBlock.Build(n, 1, n - 2, Complex.Conjugate(q), delta);
        var perm = WeightCoherenceBlock.BraComplementPermutation(n, 1, 2);
        int d = perm.Length;
        double res = 0;
        for (int t = 0; t < d; t++)
            for (int u = 0; u < d; u++)
            {
                Complex expected = -Complex.Conjugate(l12[t, u]) - (t == u ? new Complex(2.0 * n, 0) : Complex.Zero);
                res = Math.Max(res, (lpartner[perm[t], perm[u]] - expected).Magnitude);
            }
        return res;
    }

    [Fact]
    public void Zz_IsEvenUnderTheGlobalBitFlip()
    {
        // zz(c̄) = zz(c): each Z flips sign under the global bit-flip, the ZZ product is unchanged. This is the
        // structural reason the cross-fold survives the Δ·ZZ anisotropy.
        for (int n = 2; n <= 8; n++)
        {
            int full = (1 << n) - 1;
            for (int c = 0; c < (1 << n); c++)
                Assert.Equal(WeightCoherenceBlock.Zz(n, c), WeightCoherenceBlock.Zz(n, full ^ c));
        }
    }

    [Fact]
    public void Zz_KnownValues()
    {
        Assert.Equal(2, WeightCoherenceBlock.Zz(3, 0b000));   // all aligned: 2 bonds, both +1
        Assert.Equal(2, WeightCoherenceBlock.Zz(3, 0b111));   // bit-flip image of 000, same value
        Assert.Equal(-2, WeightCoherenceBlock.Zz(3, 0b010));  // 0-1 differ (−1), 1-2 differ (−1)
        Assert.Equal(0, WeightCoherenceBlock.Zz(3, 0b100));   // 0-1 equal (+1), 1-2 differ (−1)
    }

    [Theory]
    [InlineData(4)]
    [InlineData(6)]
    [InlineData(7)]
    public void DeltaZero_Overload_ReproducesThePureXyBlock(int n)
    {
        // The (q,Δ) overload at Δ=0 must be bit-identical to the legacy XY build (the delegation contract).
        var q = new Complex(1.3, -0.4);
        var xy = WeightCoherenceBlock.Build(n, 1, 2, q);
        var atZero = WeightCoherenceBlock.Build(n, 1, 2, q, 0.0);
        int d = xy.GetLength(0);
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                Assert.Equal(xy[i, j], atZero[i, j]);
    }

    [Fact]
    public void DeltaTerm_IsPresentAndDiagonal()
    {
        // Δ≠0 must actually change the block (else the test below would be vacuous), and only on the diagonal
        // (the Δ·ZZ term is diagonal in the computational basis).
        int n = 5;
        var q = new Complex(1.0, 0);
        var xy = WeightCoherenceBlock.Build(n, 1, 2, q, 0.0);
        var xxz = WeightCoherenceBlock.Build(n, 1, 2, q, 0.7);
        int d = xy.GetLength(0);
        bool anyDiagDiff = false;
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                if (i == j) anyDiagDiff |= (xxz[i, j] - xy[i, j]).Magnitude > 1e-12;
                else Assert.Equal(xy[i, j], xxz[i, j]);                       // off-diagonal untouched by Δ
        Assert.True(anyDiagDiff, "Δ=0.7 left the diagonal unchanged; the ZZ frequency term is missing");
    }

    [Theory]
    [InlineData(0.0)]      // the original F89d (XY / integrable) case
    [InlineData(0.3)]
    [InlineData(0.7)]
    [InlineData(1.0)]
    [InlineData(-0.5)]
    public void CrossFold_IsExactAntiunitarySimilarity_AtEveryDelta_RealQ(double delta)
    {
        // The headline: L(1,N−2)(q̄,Δ) = −P conj(L(1,2)(q,Δ)) Pᵀ − 2N·I to machine zero for N=4..9 at every Δ.
        // The fold is integrability-INDEPENDENT: it holds for the full interacting XXZ block, not just the
        // free-fermion XY one. (The diabolics themselves DIE under Δ; the pairing structure does not.)
        foreach (int n in SweepN)
            foreach (var qRe in new[] { 1.0, 2.0, 1.1264 })
            {
                double res = CrossFoldResidual(n, new Complex(qRe, 0), delta);
                Assert.True(res < 1e-9,
                    $"cross-fold broke at N={n}, q={qRe}, Δ={delta}: residual {res:E2}");
            }
    }

    [Theory]
    [InlineData(0.5)]
    [InlineData(1.0)]
    public void CrossFold_IsExactAntiunitarySimilarity_AtEveryDelta_ComplexQ(double delta)
    {
        // The identity is the F1 antiunitary form, so it holds at COMPLEX q too (partner at the conjugate q̄).
        foreach (int n in SweepN)
            foreach (var q in new[] { new Complex(0.5, 0.3), new Complex(0.7, -0.219) })
            {
                double res = CrossFoldResidual(n, q, delta);
                Assert.True(res < 1e-9,
                    $"cross-fold broke at N={n}, q={q}, Δ={delta}: residual {res:E2}");
            }
    }

    // The shared similarity-residual probe: max over (t,u) of
    // |partner[P t, P u] − (s(source[t,u]) − shift·δ)|, with s = −conj(·) for the antiunitary legs, identity for
    // the unitary full flip. Zero ⟹ exact similarity.
    private static double SimilarityResidual(Complex[,] source, Complex[,] partner, int[] perm, bool conjugate, double shift)
    {
        int d = perm.Length;
        double res = 0;
        for (int t = 0; t < d; t++)
            for (int u = 0; u < d; u++)
            {
                Complex s = conjugate ? -Complex.Conjugate(source[t, u]) : source[t, u];
                Complex expected = s - (t == u ? new Complex(shift, 0) : Complex.Zero);
                res = Math.Max(res, (partner[perm[t], perm[u]] - expected).Magnitude);
            }
        return res;
    }

    // Every (wKet, wBra) with both legs present and source dim ≤ cap, at the given N (a representative sweep).
    private static System.Collections.Generic.IEnumerable<(int wKet, int wBra)> Weights(int n, int dimCap = 600)
    {
        for (int wk = 1; wk < n; wk++)
            for (int wb = 1; wb < n; wb++)
                if (Binom(n, wk) * Binom(n, wb) <= dimCap)
                    yield return (wk, wb);
    }

    private static int Binom(int n, int k)
    {
        long r = 1;
        for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
        return (int)r;
    }

    [Theory]
    [InlineData(5, 0.0)]
    [InlineData(6, 0.6)]
    [InlineData(7, 1.0)]
    public void BraLeg_IsExactAntiunitarySimilarity_AtAllKetWeights(int n, double delta)
    {
        // F89d GENERALIZED past wKet=1: the bra-complement leg P (flips the bra index, right-mult ρ·F = the spine
        // V₄ element R, a factor of the F1 palindrome Π = R·D) is the EXACT antiunitary similarity
        // L(wKet,N−wBra)(q̄,Δ) = −P·conj(L(wKet,wBra)(q,Δ))·Pᵀ − 2N·I at EVERY ket weight, not just wKet=1.
        var q = new Complex(1.3, -0.2);
        foreach (var (wk, wb) in Weights(n))
        {
            var src = WeightCoherenceBlock.Build(n, wk, wb, q, delta);
            var partner = WeightCoherenceBlock.Build(n, wk, n - wb, Complex.Conjugate(q), delta);
            var perm = WeightCoherenceBlock.BraComplementPermutation(n, wk, wb);
            double res = SimilarityResidual(src, partner, perm, conjugate: true, shift: 2.0 * n);
            Assert.True(res < 1e-9, $"bra-leg broke at N={n}, (wKet,wBra)=({wk},{wb}), Δ={delta}: residual {res:E2}");
        }
    }

    [Theory]
    [InlineData(5, 0.0)]
    [InlineData(6, 0.6)]
    [InlineData(7, 1.0)]
    public void KetLeg_IsExactAntiunitarySimilarity_AtAllBraWeights(int n, double delta)
    {
        // The NEW ket-leg (the mirror of F89d on the ket index): Q (flips the ket, left-mult F·ρ = the spine V₄
        // element 𝓕R = Π²·R) gives L(N−wKet,wBra)(q̄,Δ) = −Q·conj(L(wKet,wBra)(q,Δ))·Qᵀ − 2N·I, the same
        // antiunitary form and the same −2N reflection (n_diff(ā,b) = N − n_diff(a,b) flips for the ket leg too).
        var q = new Complex(1.3, -0.2);
        foreach (var (wk, wb) in Weights(n))
        {
            var src = WeightCoherenceBlock.Build(n, wk, wb, q, delta);
            var partner = WeightCoherenceBlock.Build(n, n - wk, wb, Complex.Conjugate(q), delta);
            var perm = WeightCoherenceBlock.KetComplementPermutation(n, wk, wb);
            double res = SimilarityResidual(src, partner, perm, conjugate: true, shift: 2.0 * n);
            Assert.True(res < 1e-9, $"ket-leg broke at N={n}, (wKet,wBra)=({wk},{wb}), Δ={delta}: residual {res:E2}");
        }
    }

    [Theory]
    [InlineData(5, 0.0)]
    [InlineData(6, 0.6)]
    [InlineData(7, 1.0)]
    public void FullFlip_IsUnitarySpinFlipSimilarity_SameQ(int n, double delta)
    {
        // The full complement QP = Q∘P (the global spin-flip X^⊗N = the spine V₄ element 𝓕 = Π²): a UNITARY plain
        // similarity at the SAME q with NO conjugation and NO shift, because complementing BOTH indices leaves the
        // XOR (hence n_diff and zz) fixed: L(N−wKet,N−wBra)(q,Δ) = (QP)·L(wKet,wBra)(q,Δ)·(QP)ᵀ. This is the
        // block-resolved face of the already-typed XGlobalChargeConjugationPairing (Π²); cross-checked below.
        var q = new Complex(1.3, -0.2);
        foreach (var (wk, wb) in Weights(n))
        {
            // QP via composition: flip ket first (→ (N−wk, wb)), then flip bra (→ (N−wk, N−wb)).
            var qPerm = WeightCoherenceBlock.KetComplementPermutation(n, wk, wb);
            var pPerm = WeightCoherenceBlock.BraComplementPermutation(n, n - wk, wb);
            var full = new int[qPerm.Length];
            for (int t = 0; t < full.Length; t++) full[t] = pPerm[qPerm[t]];

            // the full-flip partner weights match the typed XGlobalChargeConjugationPairing.PairSector.
            var (pairKet, pairBra) = XGlobalChargeConjugationPairing.PairSector(n, wk, wb);
            Assert.Equal((n - wk, n - wb), (pairKet, pairBra));

            var src = WeightCoherenceBlock.Build(n, wk, wb, q, delta);
            var partner = WeightCoherenceBlock.Build(n, pairKet, pairBra, q, delta);   // SAME q, no conjugation
            double res = SimilarityResidual(src, partner, full, conjugate: false, shift: 0.0);
            Assert.True(res < 1e-9, $"full-flip (spin-flip) broke at N={n}, (wKet,wBra)=({wk},{wb}), Δ={delta}: residual {res:E2}");
        }
    }

    [Theory]
    [InlineData(5, 1, 2, 0.0)]
    [InlineData(6, 2, 3, 0.6)]
    [InlineData(6, 3, 4, 1.0)]
    public void Field_Overload_MatchesTheLongitudinalFieldHelper(int n, int wKet, int wBra, double delta)
    {
        // The promoted field overload Build(n,wKet,wBra,q,Δ,w) must equal the trusted test helper that adds the
        // longitudinal Z-field Σ_k w_k Z_k onto the (q,Δ) block — the same diagonal frequency −i·q·(fe(ket)−fe(bra)).
        // (The helper is what the filling-threshold Diagnostics harness needs in the production builder.)
        var q = new Complex(1.3, -0.4);
        double[] w = Enumerable.Range(0, n).Select(k => 0.2 * k - 0.5).ToArray();
        var expected = WithLongitudinalField(WeightCoherenceBlock.Build(n, wKet, wBra, q, delta), n, wKet, wBra, q, w);
        var actual = WeightCoherenceBlock.Build(n, wKet, wBra, q, delta, w);
        int d = expected.GetLength(0);
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                Assert.True((expected[i, j] - actual[i, j]).Magnitude < 1e-12,
                    $"field overload differs from the helper at [{i},{j}] for N={n} ({wKet},{wBra}) Δ={delta}");
    }

    [Theory]
    [InlineData(5, 1, 2, 0.5)]
    [InlineData(6, 3, 4, 0.0)]
    public void Field_Null_ReproducesTheNoFieldBuild(int n, int wKet, int wBra, double delta)
    {
        // w=null must be a no-op (the field overload reduces to the plain (q,Δ) block).
        var q = new Complex(0.9, 0.3);
        var plain = WeightCoherenceBlock.Build(n, wKet, wBra, q, delta);
        var nulled = WeightCoherenceBlock.Build(n, wKet, wBra, q, delta, null);
        int d = plain.GetLength(0);
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                Assert.Equal(plain[i, j], nulled[i, j]);
    }

    [Fact]
    public void Field_LeavesTheRealAtRateUntouched()
    {
        // The field is diagonal and Hermitian ⟹ a pure IMAGINARY frequency shift; the real Absorption-Theorem
        // rate Re λ_diag = −2·n_diff must be unchanged (the field never touches the real part on the diagonal).
        int n = 6;
        var q = new Complex(1.1, 0);
        double[] w = { 0.4, -0.3, 0.6, 0.2, -0.5, 0.1 };
        var bare = WeightCoherenceBlock.Build(n, 2, 3, q, 0.0);
        var fielded = WeightCoherenceBlock.Build(n, 2, 3, q, 0.0, w);
        int d = bare.GetLength(0);
        for (int i = 0; i < d; i++)
            Assert.Equal(bare[i, i].Real, fielded[i, i].Real, 12);   // real diagonal (AT rate) untouched
    }

    [Fact]
    public void LongitudinalZField_BreaksTheFold()
    {
        // The complementary control (so the survival result is not vacuous): a bit-flip-ODD perturbation breaks
        // the fold. A longitudinal Z-field Σ_k w_k Z_k has fieldEnergy(b̄) = −fieldEnergy(b) (each z_k flips), so
        // the bra-complement does NOT preserve its diagonal contribution, so the residual is O(1), not machine zero.
        int n = 6;
        var q = new Complex(1.3, 0);
        double[] w = { 0.4, -0.3, 0.6, 0.2, -0.5, 0.1 };
        var l12 = WithLongitudinalField(WeightCoherenceBlock.Build(n, 1, 2, q, 0.0), n, 1, 2, q, w);
        var lpartner = WithLongitudinalField(WeightCoherenceBlock.Build(n, 1, n - 2, Complex.Conjugate(q), 0.0), n, 1, n - 2, Complex.Conjugate(q), w);
        var perm = WeightCoherenceBlock.BraComplementPermutation(n, 1, 2);
        int d = perm.Length;
        double res = 0;
        for (int t = 0; t < d; t++)
            for (int u = 0; u < d; u++)
            {
                Complex expected = -Complex.Conjugate(l12[t, u]) - (t == u ? new Complex(2.0 * n, 0) : Complex.Zero);
                res = Math.Max(res, (lpartner[perm[t], perm[u]] - expected).Magnitude);
            }
        Assert.True(res > 1.0, $"a longitudinal Z-field should break the cross-fold, but residual was only {res:E2}");
    }

    // Add a longitudinal Z-field Σ_k w_k Z_k (diagonal, frequency −i·q·(fe(ket) − fe(bra)), fe(c) = Σ_k w_k·z_k,
    // z_k = −1 if site k excited else +1) on top of a (wKet,wBra) block. Used only as the negative control.
    private static Complex[,] WithLongitudinalField(Complex[,] block, int n, int wKet, int wBra, Complex q, double[] w)
    {
        var kets = WeightCoherenceBlock.Configs(n, wKet);
        var bras = WeightCoherenceBlock.Configs(n, wBra);
        int col = 0;
        foreach (var kc in kets)
            foreach (var bc in bras)
            {
                double fe = FieldEnergy(n, w, kc) - FieldEnergy(n, w, bc);
                block[col, col] += (-Complex.ImaginaryOne) * q * fe;
                col++;
            }
        return block;
    }

    private static double FieldEnergy(int n, double[] w, int c)
    {
        double e = 0;
        for (int k = 0; k < n; k++) e += w[k] * (((c >> k) & 1) == 1 ? -1.0 : 1.0);
        return e;
    }
}
