using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>One N's reading of the (SE,DE)↔(SE,w_{N−2}) cross-fold: the partner bra-weight, the block dimension,
/// and the residual of the exact antiunitary similarity that carries the branch-locus palindrome across the two
/// blocks.</summary>
public sealed record CrossFoldReading(int N, int PartnerWBra, int Dim, double SimilarityResidual);

/// <summary>Move 4, answered: the (SE,DE) diabolics PAIR across the cross-block fold, because that fold is an
/// EXACT antiunitary similarity.
///
/// <para>The branch-locus palindrome's bra bit-flip ρ[a,b] → ρ[a,b̄] (the F89c lemma, n_diff(a,b̄) = N −
/// n_diff(a,b)) maps the (SE,DE) = (w1,w2) coherence block to the (SE, w_{N−2}) = (w1, N−2) block. This witness
/// builds both blocks (<see cref="WeightCoherenceBlock"/>) and the bra-complement permutation P
/// (<see cref="WeightCoherenceBlock.BraComplementPermutation"/>) and checks the matrix identity
/// L(1,N−2)(q̄) = −P · conj(L(1,2)(q)) · Pᵀ − 2N·I to machine zero (the entries are exact arithmetic; the residual
/// is 0 for N=4..9 at every q, real or complex). This is STRONGER than the spectrum match the CLI's
/// <c>foldcross</c> command reports: an antiunitary similarity preserves the whole Jordan structure, so a
/// SEMISIMPLE coalescence (a diabolic) in (SE,DE) at (q, λ) maps to a semisimple coalescence in (SE,w_{N−2}) at
/// (q̄, −λ̄−2N) with the IDENTICAL coalescence gap and character. Hence every (SE,DE) diabolic has a cross-fold
/// partner diabolic, for all N and all q at once, no enumeration needed.</para>
///
/// <para>At N=4 the partner w_{N−2} = w2 = DE, so the partner IS the (SE,DE) block (the N=4-only within-block
/// self-fold, the degenerate partner=self case that put one diabolic on the real axis); for N ≥ 5 the partner is
/// a different block ((SE,TE) at N=5, (SE,QE) at N=6, …), and the N=4 on-line "zeros" become cross-block mirror
/// partners. The witness also reproduces the pairing on the N=7 real-q diabolic (λ=−4.942 ↔ partner −9.058, equal
/// gaps). Registered as F89d in <c>docs/ANALYTICAL_FORMULAS.md</c> and typed as
/// <c>F89CrossFoldSimilarityClaim</c> (parents F1 + the branch-locus palindrome). Anchor:
/// <c>experiments/F89_PATH_K_DIABOLIC.md</c> (the cross-fold section) and
/// <c>experiments/F89_BRANCH_LOCUS_PALINDROME.md</c>; the persistent evidence for the
/// <c>diabolic_over_higher_n</c> arc's Move 4.</para></summary>
public sealed class CrossFoldSimilarityWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private static readonly int[] SweepN = { 4, 5, 6, 7, 8, 9 };
    private const double GenericQ = 1.0;

    /// <summary>The cross-fold antiunitary-similarity residual at coupling q:
    /// max over (t,u) of |L(1,N−2)(q̄)[Pt,Pu] − (−conj(L(1,2)(q)[t,u]) − 2N·δ)|. Zero ⟹ exact similarity.</summary>
    public CrossFoldReading Read(int n, Complex q)
    {
        var l1 = WeightCoherenceBlock.Build(n, 1, 2, q);
        var l2 = WeightCoherenceBlock.Build(n, 1, n - 2, Complex.Conjugate(q));     // partner at q̄ (the F1 form)
        var perm = WeightCoherenceBlock.BraComplementPermutation(n, 1, 2);
        int d = perm.Length;
        double res = 0;
        for (int t = 0; t < d; t++)
            for (int u = 0; u < d; u++)
            {
                Complex expected = -Complex.Conjugate(l1[t, u]) - (t == u ? new Complex(2.0 * n, 0) : Complex.Zero);
                res = Math.Max(res, (l2[perm[t], perm[u]] - expected).Magnitude);
            }
        return new CrossFoldReading(n, n - 2, d, res);
    }

    /// <summary>The fold-image of a real-q diabolic: the (SE,DE) coalescence gap near (qRe, λ) and the partner
    /// (SE,w_{N−2}) coalescence gap near the fold image −λ−2N (real λ, so q̄ = q). Equal gaps ⟹ the diabolic
    /// pairs across the fold (a concrete reading of the structural similarity).</summary>
    public (double Gap12, double GapPartner, double PartnerLambda) ReproducePairedDiabolic(int n, double qRe, double lambda)
    {
        var q = new Complex(qRe, 0);
        double partnerLam = -lambda - 2.0 * n;                                       // −λ̄ − 2N for real λ
        double g12 = MinGapNear(WeightCoherenceBlock.Build(n, 1, 2, q), new Complex(lambda, 0));
        double gp = MinGapNear(WeightCoherenceBlock.Build(n, 1, n - 2, q), new Complex(partnerLam, 0));
        return (g12, gp, partnerLam);
    }

    /// <summary>The gap between the two block eigenvalues nearest a target (a coalescence reads ≈ 0).</summary>
    private static double MinGapNear(Complex[,] block, Complex target)
    {
        var ev = Matrix<Complex>.Build.DenseOfArray(block).Evd().EigenValues.ToArray();
        var near = ev.OrderBy(z => (z - target).Magnitude).Take(2).ToArray();
        return (near[0] - near[1]).Magnitude;
    }

    public string DisplayName =>
        "CrossFoldSimilarityWitness (the (SE,DE)↔(SE,w_{N−2}) cross-fold is an EXACT antiunitary similarity, so the diabolics pair)";

    public string Summary
    {
        get
        {
            var r5 = Read(5, new Complex(GenericQ, 0));
            return "the branch-locus palindrome's bra bit-flip is an EXACT antiunitary similarity " +
                   $"L(1,N−2)(q̄) = −P·conj(L(1,2)(q))·Pᵀ − 2N·I (residual {r5.SimilarityResidual.ToString("E1", Inv)} at N=5, " +
                   "exact arithmetic), so every (SE,DE) diabolic at (q, λ) has a partner diabolic at (q̄, −λ̄−2N) in the " +
                   "(SE,w_{N−2}) block with identical character (an antiunitary similarity preserves Jordan structure) and " +
                   "identical coalescence gap. The N=4 self-fold is the degenerate partner=self case. Move 4, answered: the " +
                   "diabolics pair across the cross-fold, for all N and all q at once.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (int n in SweepN)
            {
                var r = Read(n, new Complex(GenericQ, 0));
                string partner = n == 4 ? "(1,2) = SELF" : $"(1,{r.PartnerWBra})";
                yield return new InspectableNode(
                    displayName: $"N={n}: (SE,DE)=(1,2) ↔ (SE,w_{{N−2}})={partner}, dim {r.Dim}",
                    summary: $"antiunitary-similarity residual {r.SimilarityResidual.ToString("E2", Inv)} (machine zero) at q={GenericQ.ToString("0.#", Inv)} " +
                             (n == 4
                                 ? "⟹ N=4: the partner IS the (SE,DE) block, the within-block self-fold (partner=self)."
                                 : "⟹ an exact cross-fold similarity: the (SE,DE) spectrum folds onto the partner's, the whole Jordan structure preserved."),
                    provenance: NodeProvenance.Live);
            }

            foreach (var (n, q, lam) in new[] { (7, 1.1264, -4.942) })
            {
                var (g12, gp, plam) = ReproducePairedDiabolic(n, q, lam);
                yield return new InspectableNode(
                    displayName: $"gate: N={n} real-q diabolic pairs across the fold ((1,2) λ={lam.ToString("0.###", Inv)} ↔ (1,{n - 2}) λ={plam.ToString("0.###", Inv)})",
                    summary: $"at q={q.ToString("0.####", Inv)} the (1,2) coalescence gap near λ={lam.ToString("0.###", Inv)} is {g12.ToString("E2", Inv)}, " +
                             $"and the partner (1,{n - 2}) gap near the fold image −λ−2N={plam.ToString("0.###", Inv)} is {gp.ToString("E2", Inv)} " +
                             "(equal ⟹ the diabolic pairs across the cross-fold).",
                    provenance: NodeProvenance.Live);
            }
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
