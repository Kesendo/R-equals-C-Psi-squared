using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>One N's reading of the (SE,DE)↔(SE,w_{N−2}) cross-fold at anisotropy Δ: the partner bra-weight, the
/// block dimension, and the residual of the exact antiunitary similarity that carries the branch-locus palindrome
/// across the two blocks. The residual is machine-zero at EVERY Δ (the fold is integrability-independent).</summary>
public sealed record CrossFoldReading(int N, int PartnerWBra, int Dim, double SimilarityResidual, double Delta);

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
/// <para>The fold is integrability-INDEPENDENT: the identity holds for the FULL interacting XXZ block at EVERY
/// anisotropy Δ (<see cref="WeightCoherenceBlock.Build(int,int,int,System.Numerics.Complex,double)"/>), because the
/// Δ·ZZ term is EVEN under the global bit-flip (zz(b̄) = zz(b)), so the bra-complement carries it cleanly. The
/// diabolics themselves DIE under Δ (integrability-protected, the arc's Move 2), but the pairing structure does
/// not: a diabolic and its cross-fold partner turn defective in lockstep. The discriminant is bit-flip PARITY: a
/// bit-flip-ODD perturbation breaks the fold; a longitudinal Z-field Σ_k w_k Z_k has fe(b̄) = −fe(b), so its
/// residual is O(1), not machine zero (<see cref="ReadFieldControlResidual"/>, the complementary control). The fold
/// is therefore a structural/algebraic property of the Liouvillian, not a free-fermion artifact.</para>
///
/// <para>At N=4 the partner w_{N−2} = w2 = DE, so the partner IS the (SE,DE) block (the N=4-only within-block
/// self-fold, the degenerate partner=self case that put one diabolic on the real axis); for N ≥ 5 the partner is
/// a different block ((SE,TE) at N=5, (SE,QE) at N=6, …), and the N=4 on-line "zeros" become cross-block mirror
/// partners. The witness also reproduces the pairing on the N=7 real-q diabolic (λ=−4.942 ↔ partner −9.058, equal
/// gaps).</para>
///
/// <para>The cross-fold is one leg of a Klein four-group of bit-flip similarities on the coherence-block lattice,
/// general in BOTH weights (not just wKet=1): the bra-complement P (flip the bra, <see cref="BraLegResidual"/>)
/// and the ket-complement Q (flip the ket, <see cref="KetLegResidual"/>) are the two exact ANTIUNITARY legs (each
/// with the −2N reflection), and their product is the UNITARY global spin-flip QP = X^⊗N
/// (<see cref="FullFlipResidual"/>, same q, no shift). This is not a new group: it is the existing spine
/// V₄ = {I, F⊗F, I⊗F, F⊗I} ⊂ D₄ (PROOF_PI_FACTORS_AS_R_TIMES_D / F118 MirrorGroupD4Claim), here block-resolved and
/// q-parameterized. The dock onto the F1 palindrome trunk is exact: the bra leg P = right-mult ρ·F = the spine R
/// is a FACTOR of Π (Π = R·D, D = transpose), the full flip QP = Π² = the typed
/// <c>XGlobalChargeConjugationPairing</c>, and the −2N shift is the block image of R·L_diss·R = −L_diss − 2σ.
/// Naming bridge: F89 names the legs by the flipped index (bra/ket); the D₄ proof docs name them by the
/// multiplication side, calling ρ·F the "ket reflection" (the opposite word for the same operator).</para>
///
/// <para>Registered as F89d in <c>docs/ANALYTICAL_FORMULAS.md</c> and typed as
/// <c>F89CrossFoldSimilarityClaim</c> (parents F1 + the branch-locus palindrome). Anchor:
/// <c>experiments/F89_PATH_K_DIABOLIC.md</c> (the cross-fold section) and
/// <c>experiments/F89_BRANCH_LOCUS_PALINDROME.md</c>; the persistent evidence for the
/// <c>diabolic_over_higher_n</c> arc's Move 4.</para></summary>
public sealed class CrossFoldSimilarityWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private static readonly int[] SweepN = { 4, 5, 6, 7, 8, 9 };
    private const double GenericQ = 1.0;

    /// <summary>The cross-fold antiunitary-similarity residual of the pure-XY (Δ=0) block at coupling q.</summary>
    public CrossFoldReading Read(int n, Complex q) => Read(n, q, 0.0);

    /// <summary>The cross-fold antiunitary-similarity residual of the XXZ block at coupling q and anisotropy Δ:
    /// max over (t,u) of |L(1,N−2)(q̄,Δ)[Pt,Pu] − (−conj(L(1,2)(q,Δ)[t,u]) − 2N·δ)|. Zero ⟹ exact similarity.
    /// Machine-zero at EVERY Δ: the Δ·ZZ term is even under the global bit-flip (zz(b̄) = zz(b)), so the bra-
    /// complement carries it cleanly; the fold is integrability-independent (it survives the anisotropy that
    /// kills the diabolics themselves).</summary>
    public CrossFoldReading Read(int n, Complex q, double delta)
    {
        var l1 = WeightCoherenceBlock.Build(n, 1, 2, q, delta);
        var l2 = WeightCoherenceBlock.Build(n, 1, n - 2, Complex.Conjugate(q), delta);   // partner at q̄ (the F1 form)
        var perm = WeightCoherenceBlock.BraComplementPermutation(n, 1, 2);
        int d = perm.Length;
        double res = 0;
        for (int t = 0; t < d; t++)
            for (int u = 0; u < d; u++)
            {
                Complex expected = -Complex.Conjugate(l1[t, u]) - (t == u ? new Complex(2.0 * n, 0) : Complex.Zero);
                res = Math.Max(res, (l2[perm[t], perm[u]] - expected).Magnitude);
            }
        return new CrossFoldReading(n, n - 2, d, res, delta);
    }

    /// <summary>The cross-fold residual with a longitudinal Z-field Σ_k w_k Z_k added to BOTH blocks (the diagonal
    /// frequency −i·q·(fe(ket) − fe(bra)), fe(c) = Σ_k w_k·z_k, z_k = −1 if site k excited else +1). The
    /// COMPLEMENTARY control: a field is ODD under the global bit-flip (fe(b̄) = −fe(b)), so the bra-complement
    /// flips its sign and the residual is O(1), NOT machine zero, pinning the discriminant as bit-flip parity
    /// (even ZZ survives, odd field breaks).</summary>
    public double ReadFieldControlResidual(int n, Complex q, double[] field)
    {
        var l1 = WithLongitudinalField(WeightCoherenceBlock.Build(n, 1, 2, q), n, 1, 2, q, field);
        var l2 = WithLongitudinalField(WeightCoherenceBlock.Build(n, 1, n - 2, Complex.Conjugate(q)), n, 1, n - 2, Complex.Conjugate(q), field);
        var perm = WeightCoherenceBlock.BraComplementPermutation(n, 1, 2);
        int d = perm.Length;
        double res = 0;
        for (int t = 0; t < d; t++)
            for (int u = 0; u < d; u++)
            {
                Complex expected = -Complex.Conjugate(l1[t, u]) - (t == u ? new Complex(2.0 * n, 0) : Complex.Zero);
                res = Math.Max(res, (l2[perm[t], perm[u]] - expected).Magnitude);
            }
        return res;
    }

    /// <summary>The bra-leg residual at GENERAL (wKet, wBra): |L(wKet,N−wBra)(q̄,Δ)[Pt,Pu] − (−conj(L(wKet,wBra)
    /// (q,Δ)[t,u]) − 2N·δ)|. F89d generalized past wKet=1; P = the bra-complement (right-mult ρ·F = the spine V₄
    /// element R, a factor of Π = R·D). Machine zero at every ket weight.</summary>
    public double BraLegResidual(int n, int wKet, int wBra, Complex q, double delta)
        => LegResidual(WeightCoherenceBlock.Build(n, wKet, wBra, q, delta),
                       WeightCoherenceBlock.Build(n, wKet, n - wBra, Complex.Conjugate(q), delta),
                       WeightCoherenceBlock.BraComplementPermutation(n, wKet, wBra), conjugate: true, shift: 2.0 * n);

    /// <summary>The KET-leg residual (the mirror of F89d on the ket index): |L(N−wKet,wBra)(q̄,Δ)[Qt,Qu] −
    /// (−conj(L(wKet,wBra)(q,Δ)[t,u]) − 2N·δ)|. Q = the ket-complement (left-mult F·ρ = the spine V₄ element
    /// 𝓕R = Π²·R). Same antiunitary form, same −2N reflection. Machine zero at every bra weight.</summary>
    public double KetLegResidual(int n, int wKet, int wBra, Complex q, double delta)
        => LegResidual(WeightCoherenceBlock.Build(n, wKet, wBra, q, delta),
                       WeightCoherenceBlock.Build(n, n - wKet, wBra, Complex.Conjugate(q), delta),
                       WeightCoherenceBlock.KetComplementPermutation(n, wKet, wBra), conjugate: true, shift: 2.0 * n);

    /// <summary>The full-flip residual (the global spin-flip QP = X^⊗N = Π²): |L(N−wKet,N−wBra)(q,Δ)[QPt,QPu] −
    /// L(wKet,wBra)(q,Δ)[t,u]|. A UNITARY plain similarity at the SAME q, no conjugation, no shift (complementing
    /// both indices leaves n_diff and zz fixed). The block-resolved face of XGlobalChargeConjugationPairing.</summary>
    public double FullFlipResidual(int n, int wKet, int wBra, Complex q, double delta)
    {
        var qPerm = WeightCoherenceBlock.KetComplementPermutation(n, wKet, wBra);
        var pPerm = WeightCoherenceBlock.BraComplementPermutation(n, n - wKet, wBra);
        var full = new int[qPerm.Length];
        for (int t = 0; t < full.Length; t++) full[t] = pPerm[qPerm[t]];        // QP = P∘Q
        return LegResidual(WeightCoherenceBlock.Build(n, wKet, wBra, q, delta),
                           WeightCoherenceBlock.Build(n, n - wKet, n - wBra, q, delta), full, conjugate: false, shift: 0.0);
    }

    private static double LegResidual(Complex[,] source, Complex[,] partner, int[] perm, bool conjugate, double shift)
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

    private static Complex[,] WithLongitudinalField(Complex[,] block, int n, int wKet, int wBra, Complex q, double[] w)
    {
        var kets = WeightCoherenceBlock.Configs(n, wKet);
        var bras = WeightCoherenceBlock.Configs(n, wBra);
        int col = 0;
        foreach (var kc in kets)
            foreach (var bc in bras)
            {
                double fe = FieldEnergy(n, w, kc) - FieldEnergy(n, w, bc);
                block[col, col] += (-Complex.ImaginaryOne) * q * fe;     // same (ket-outer, bra-inner) order as Build
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
            var r5d = Read(5, new Complex(GenericQ, 0), 0.7);
            return "the branch-locus palindrome's bra bit-flip is an EXACT antiunitary similarity " +
                   $"L(1,N−2)(q̄,Δ) = −P·conj(L(1,2)(q,Δ))·Pᵀ − 2N·I (residual {r5.SimilarityResidual.ToString("E1", Inv)} at N=5 Δ=0, " +
                   $"{r5d.SimilarityResidual.ToString("E1", Inv)} at N=5 Δ=0.7, exact arithmetic), so every (SE,DE) diabolic at " +
                   "(q, λ) has a partner diabolic at (q̄, −λ̄−2N) in the (SE,w_{N−2}) block with identical character (an " +
                   "antiunitary similarity preserves Jordan structure) and identical coalescence gap. The N=4 self-fold is the " +
                   "degenerate partner=self case. The fold is integrability-INDEPENDENT: it survives at EVERY Δ (the Δ·ZZ term " +
                   "is even under the global bit-flip, zz(b̄)=zz(b)), so it holds for the full interacting XXZ block even though " +
                   "the diabolics themselves die under Δ; the discriminant is bit-flip parity (a longitudinal Z-field, odd, " +
                   "breaks it). Move 4, answered: the diabolics pair across the cross-fold, for all N, all q, all Δ at once. " +
                   "The fold holds at EVERY ket weight (not just wKet=1), and has a mirror KET leg (flip the ket index): both " +
                   "legs are exact antiunitary similarities (−2N), their product the unitary global spin-flip. These are the " +
                   "existing spine V₄ ⊂ D₄ block-resolved: the bra leg P = ρ·F is a factor of the F1 palindrome Π = R·D, " +
                   "QP = Π² = XGlobalChargeConjugationPairing. So F89d docks onto F1.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (int n in SweepN)
            {
                var r = Read(n, new Complex(GenericQ, 0));
                var rd = Read(n, new Complex(GenericQ, 0), 0.7);                 // the XXZ (Δ≠0) reading
                string partner = n == 4 ? "(1,2) = SELF" : $"(1,{r.PartnerWBra})";
                yield return new InspectableNode(
                    displayName: $"N={n}: (SE,DE)=(1,2) ↔ (SE,w_{{N−2}})={partner}, dim {r.Dim}",
                    summary: $"antiunitary-similarity residual {r.SimilarityResidual.ToString("E2", Inv)} at Δ=0 and " +
                             $"{rd.SimilarityResidual.ToString("E2", Inv)} at Δ=0.7 (both machine zero) at q={GenericQ.ToString("0.#", Inv)} " +
                             (n == 4
                                 ? "⟹ N=4: the partner IS the (SE,DE) block, the within-block self-fold (partner=self)."
                                 : "⟹ an exact cross-fold similarity at every Δ: the (SE,DE) spectrum folds onto the partner's, the whole Jordan structure preserved."),
                    provenance: NodeProvenance.Live);
            }

            // The Klein four-group of bit-flip similarities (general weight): F89d (the bra leg) is one of two
            // antiunitary legs; the other is the ket leg; their product is the unitary global spin-flip.
            var qg = new Complex(1.3, -0.2);
            double braGen = BraLegResidual(6, 2, 3, qg, 0.6);                  // bra leg at wKet=2 (past wKet=1)
            double ketGen = KetLegResidual(6, 2, 3, qg, 0.6);                  // the NEW ket leg
            double fullGen = FullFlipResidual(6, 2, 3, qg, 0.6);              // the unitary spin-flip
            yield return new InspectableNode(
                displayName: "the Klein four-group of bit-flip similarities (N=6, (wKet,wBra)=(2,3), q=1.3−0.2i, Δ=0.6)",
                summary: $"bra-leg P (flips bra, ρ·F): (2,3)→(2,3) residual {braGen.ToString("E2", Inv)}; ket-leg Q (flips ket, F·ρ): " +
                         $"(2,3)→(4,3) residual {ketGen.ToString("E2", Inv)}; full QP (spin-flip X^⊗N): (2,3)→(4,3) residual {fullGen.ToString("E2", Inv)}. " +
                         "Both legs are exact antiunitary similarities (−2N reflection); F89d is the wKet=1 corner of the bra leg, " +
                         "the ket leg is its mirror, and QP=P∘Q is the unitary global spin-flip (same q, no shift, n_diff preserved).",
                provenance: NodeProvenance.Live);

            // The dock onto the F1 palindrome trunk: the legs are block-restrictions of the spine V₄ ⊂ D₄.
            yield return new InspectableNode(
                displayName: "dock: the legs are block-restrictions of the spine V₄ ⊂ D₄, factoring the F1 palindrome Π",
                summary: "the three operators ARE the existing spine V₄ = {I, F⊗F, I⊗F, F⊗I} (F=X^⊗N), block-resolved and " +
                         "q-parameterized: the bra leg P = right-mult ρ·F = the spine R, a FACTOR of Π (Π = R·D, D = transpose); " +
                         "the full flip QP = F⊗F = Π² = the typed XGlobalChargeConjugationPairing; the ket leg Q = left-mult F·ρ = " +
                         "𝓕R = Π²·R. So F89d is not a one-off: it is the F1 palindrome's bra leg, its −2N shift the block image of " +
                         "R·L_diss·R = −L_diss − 2σ (PROOF_PI_FACTORS_AS_R_TIMES_D / F118). Naming bridge: F89 names by the flipped " +
                         "index (bra/ket); the D₄ docs name by the multiplication side, calling ρ·F the 'ket reflection' (opposite word).",
                provenance: NodeProvenance.Stored);

            // The bit-flip-parity discriminant: ZZ anisotropy (even) survives, a longitudinal Z-field (odd) breaks.
            double[] zField = { 0.4, -0.3, 0.6, 0.2, -0.5, 0.1 };              // N=6 random per-site field
            double fieldRes = ReadFieldControlResidual(6, new Complex(1.3, 0), zField);
            yield return new InspectableNode(
                displayName: "the discriminant is bit-flip parity (even ZZ survives, odd Z-field breaks)",
                summary: "the cross-fold survives any bit-flip-EVEN bond term (the Δ·ZZ anisotropy, so the FULL interacting XXZ " +
                         "chain, not just the integrable XY one: zz(b̄)=zz(b)), but BREAKS under a bit-flip-ODD term: a " +
                         $"longitudinal Z-field Σ_k w_k Z_k gives residual {fieldRes.ToString("E2", Inv)} (O(1), not machine zero) at N=6 " +
                         "because fe(b̄)=−fe(b). So the fold is structural/algebraic, NOT a free-fermion (integrability) artifact.",
                provenance: NodeProvenance.Live);

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
