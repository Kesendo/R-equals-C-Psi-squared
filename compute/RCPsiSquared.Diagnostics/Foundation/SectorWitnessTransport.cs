using System;
using System.Collections.Generic;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The MEMBER instrument of the step-3 shell census (Task 4, sectorbraid large-N program;
/// docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md §5-§7). At a diamond member cell the eigenvalue sits AT
/// the probe shift, so σ_min ~ 1e-13 and κ ~ 1e14: iterative solvers cannot read it. Instead the
/// spectator intertwiner W (Theorem B; F125) transports the (1,2) block's near-defective 2-plane up
/// the |p−q̃|=1 band, the §7 fold gauge P𝒟 carries it to the μ = −λ_A−2N blocks, and a minimized
/// 2-plane residual (Householder QR of the d×2 matrix (L−s)Q, then a stable 2×2 SVD) yields a
/// certified FROM-ABOVE upper bound on σ_min with two sparse matvecs and NO solves.
///
/// <para><b>Conventions.</b> Everything runs in the ONE pinned little-endian
/// <see cref="WeightCoherenceBlock.Configs"/> coordinate system (bit s = site s, full-block flat index
/// = ketIndex·nBra + braIndex), matching the census CSR (<see cref="WeightCoherenceSectorCsr"/>). ApplyW
/// is the element rule of <see cref="SpectatorIntertwiner.BuildW"/> (Jordan-Wigner signs; that file is
/// the source of truth) TRANSLATED into this convention: BuildW's big-endian JW string over the sites
/// ABOVE bit N−1−l becomes here the string over the occupied bits strictly BELOW bit l. The arbiter of
/// correctness is the CSR intertwining gate (ApplyW_Intertwines_TheCensusCsrPencil), not literal
/// similarity to BuildW.</para>
///
/// <para><b>The fold gauge.</b> The §7 bra fold is L(p,N−q̃)(q) = −P𝒟·L(p,q̃)(q)·𝒟Pᵀ − 2N·I at every q
/// (Δ=0), with P = <see cref="WeightCoherenceBlock.BraComplementPermutation"/> (b ↦ b̄, mapping
/// (p,q̃)→(p,N−q̃)) and 𝒟 the bipartite alternating-site-sign diagonal 𝒟|a⟩⟨b| = χ(a)χ(b)|a⟩⟨b|,
/// χ(c) = (−1)^(occupied sites on the odd sublattice). <see cref="ApplyFoldGauge"/> is v ↦ (P·𝒟)v (apply
/// 𝒟 first, then P), a real signed permutation, hence an isometry; the shifted-residual isometry
/// (FoldGauge_PreservesShiftedResidual) is the pin of the exact convention.</para>
///
/// <para><b>R-parity alternation (F7, §7).</b> R·W = (−1)^(p+q̃)·W·R, so every band step off an
/// odd-(p+q̃) block FLIPS the carried R-parity and the folds preserve it: carriedOdd(band member
/// (p,p+1)) = seedOdd XOR ((p−1) mod 2).</para></summary>
public static class SectorWitnessTransport
{
    // -----------------------------------------------------------------------------------------------
    // W as a matrix-free action (p,q̃) -> (p+1,q̃+1)
    // -----------------------------------------------------------------------------------------------

    /// <summary>W(ρ) = Σ_l c_l†ρc_l as a matrix-free action (p,w)→(p+1,w+1) in little-endian
    /// <see cref="WeightCoherenceBlock.Configs"/> coordinates (row = ketIndex·nBra + braIndex). For each
    /// input coherence |a⟩⟨b| and each site l unoccupied in BOTH a and b, adds l to bra and ket with
    /// coefficient s_l(a)·s_l(b), s_l(c) = (−1)^(occupied bits strictly below l) the little-endian
    /// Jordan-Wigner string (<see cref="SpectatorIntertwiner.BuildW"/> mirrored).</summary>
    public static Complex[] ApplyW(int n, int p, int w, Complex[] v)
    {
        if (v is null) throw new ArgumentNullException(nameof(v));
        var kets = WeightCoherenceBlock.Configs(n, p);
        var bras = WeightCoherenceBlock.Configs(n, w);
        int nBra = bras.Count;
        if (v.Length != kets.Count * nBra) throw new ArgumentException($"v length {v.Length} != dim {kets.Count * nBra}", nameof(v));

        var ketsOut = WeightCoherenceBlock.Configs(n, p + 1);
        var brasOut = WeightCoherenceBlock.Configs(n, w + 1);
        int nBraOut = brasOut.Count;
        var ketIndexOut = Index(ketsOut);
        var braIndexOut = Index(brasOut);

        var result = new Complex[ketsOut.Count * nBraOut];
        for (int col = 0; col < v.Length; col++)
        {
            Complex vc = v[col];
            if (vc == Complex.Zero) continue;
            int a = kets[col / nBra], b = bras[col % nBra];
            int occupied = a | b;
            for (int l = 0; l < n; l++)
            {
                int bit = 1 << l;
                if ((occupied & bit) != 0) continue;                       // l must be empty in both a and b
                int sign = JwSign(a, l) * JwSign(b, l);
                int row = ketIndexOut[a | bit] * nBraOut + braIndexOut[b | bit];
                result[row] += sign * vc;
            }
        }
        return result;
    }

    /// <summary>The little-endian Jordan-Wigner string sign s_l(c) = (−1)^(# occupied bits strictly
    /// below bit l). Mirror of <see cref="SpectatorIntertwiner.JwSign"/> (which counts the occupied
    /// sites before l = the bits ABOVE bit N−1−l in the big-endian PopcountStates basis).</summary>
    private static int JwSign(int c, int l)
        => (BitOperations.PopCount((uint)(c & ((1 << l) - 1))) & 1) == 0 ? 1 : -1;

    // -----------------------------------------------------------------------------------------------
    // the §7 fold gauge v -> (P·𝒟)v : (p,w) -> (p, N-w)
    // -----------------------------------------------------------------------------------------------

    /// <summary>The §7 fold gauge v ↦ (P·𝒟)v from (p,w) coordinates to (p,N−w) coordinates, a real
    /// signed permutation (hence an isometry). 𝒟 is the bipartite site-alternating sign
    /// χ(ket)·χ(bra); P = <see cref="WeightCoherenceBlock.BraComplementPermutation"/> re-indexes the
    /// bra to its complement. Composing the two realizes the left factor of the bra fold identity
    /// L(p,N−w) = −(P𝒟)·L(p,w)·(𝒟Pᵀ) − 2N·I.</summary>
    public static Complex[] ApplyFoldGauge(int n, int p, int w, Complex[] v)
    {
        if (v is null) throw new ArgumentNullException(nameof(v));
        var kets = WeightCoherenceBlock.Configs(n, p);
        var bras = WeightCoherenceBlock.Configs(n, w);
        int nBra = bras.Count;
        if (v.Length != kets.Count * nBra) throw new ArgumentException($"v length {v.Length} != dim {kets.Count * nBra}", nameof(v));
        var perm = WeightCoherenceBlock.BraComplementPermutation(n, p, w);   // (p,w) index -> (p,N-w) index
        int oddMask = OddSublatticeMask(n);

        var result = new Complex[v.Length];                                 // dim(p,N-w) = dim(p,w)
        for (int t = 0; t < v.Length; t++)
        {
            int a = kets[t / nBra], b = bras[t % nBra];
            int sign = Chi(a, oddMask) * Chi(b, oddMask);                   // 𝒟 acts first
            result[perm[t]] = sign * v[t];                                  // then P re-indexes
        }
        return result;
    }

    /// <summary>The bipartite gauge sign χ(c) = (−1)^(occupied bits on the odd sublattice) = ⟨c|D|c⟩ for
    /// D = ∏_(j odd) Z_j (D H_XY D = −H_XY on the open chain). The choice of sublattice is silent for
    /// the fold identity: swapping odd↔even multiplies 𝒟 by a block-global sign that cancels between
    /// the two 𝒟 factors.</summary>
    private static int Chi(int c, int oddMask)
        => (BitOperations.PopCount((uint)(c & oddMask)) & 1) == 0 ? 1 : -1;

    private static int OddSublatticeMask(int n)
    {
        int m = 0;
        for (int s = 1; s < n; s += 2) m |= 1 << s;
        return m;
    }

    // -----------------------------------------------------------------------------------------------
    // from-above member reading: Householder 2-plane residual [N1]
    // -----------------------------------------------------------------------------------------------

    /// <summary>An UPPER BOUND on σ_min(L−s) from the 2-plane spanned by {v1,v2}: orthonormalize the
    /// plane (2-column modified Gram-Schmidt), take two sparse shifted matvecs w_i = (L−s)q_i, then a
    /// 2-column complex Householder QR of the d×2 matrix [w1 w2] followed by a STABLE 2×2 SVD [N1]. The
    /// value equals min over unit x in the plane of ‖(L−s)x‖ = σ_min([w1 w2]) ≥ σ_min(L−s). NO Gram
    /// matrix Qᴴ(L−s)ᴴ(L−s)Q anywhere (that squares the member depth below the rounding floor); σ_min
    /// of the 2×2 upper-triangular R = [[a,b],[0,c]] is read as |a·c|/σ_max, σ_max the large SVD root,
    /// which never subtracts nearly-equal squares.</summary>
    public static double TwoPlaneResidual(WeightCoherenceSectorCsr.Csr block, Complex shift, Complex[] v1, Complex[] v2)
    {
        if (block is null) throw new ArgumentNullException(nameof(block));
        int d = block.Dim;
        if (v1 is null || v2 is null) throw new ArgumentNullException(v1 is null ? nameof(v1) : nameof(v2));
        if (v1.Length != d || v2.Length != d) throw new ArgumentException("plane vectors must have length block.Dim");

        // 2-column orthonormalization of {v1,v2} (well-conditioned invariant subspace; [N3])
        var q1 = (Complex[])v1.Clone();
        var q2 = (Complex[])v2.Clone();
        double n1 = Norm(q1);
        if (n1 == 0) return 0.0;
        Scale(q1, 1.0 / n1);
        Complex proj = Inner(q1, q2);                                       // <q1,q2>
        for (int i = 0; i < d; i++) q2[i] -= proj * q1[i];
        double n2 = Norm(q2);
        if (n2 == 0) return 0.0;                                            // degenerate plane
        Scale(q2, 1.0 / n2);

        // two sparse shifted matvecs: M = [w1 w2] = (L - s) Q
        var w1 = new Complex[d];
        var w2 = new Complex[d];
        CsrOps.MultiplyShifted(block, shift, q1, w1);
        CsrOps.MultiplyShifted(block, shift, q2, w2);

        // 2-column complex Householder QR of the d x 2 matrix [w1 w2] -> 2x2 upper triangular R
        Complex r00 = HouseholderColumn(w1, 0, d, out Complex[]? hv0, out double beta0);
        ApplyReflector(hv0, beta0, w2, 0, d);
        Complex r01 = w2[0];
        Complex r11 = HouseholderColumn(w2, 1, d, out _, out _);

        return SigmaMin2x2(r00, r01, r11);
    }

    /// <summary>Stable smallest singular value of the 2×2 upper-triangular R = [[a,b],[0,c]]. det R = a·c,
    /// so σ_max·σ_min = |a·c|; with the Frobenius sum F = |a|²+|b|²+|c|², σ_max² = (F + √(F²−4|ac|²))/2
    /// (the large root, no cancellation) and σ_min = |ac|/σ_max (a division, stable when σ_min ≪ σ_max).
    /// This is the [N1] route that resolves the quadratic member depth the Gram matrix would lose.</summary>
    private static double SigmaMin2x2(Complex a, Complex b, Complex c)
    {
        double detAbs = a.Magnitude * c.Magnitude;                         // |det R| = |a|·|c|
        double f = a.Magnitude * a.Magnitude + b.Magnitude * b.Magnitude + c.Magnitude * c.Magnitude;
        double disc = f * f - 4.0 * detAbs * detAbs;
        double sqrtDisc = Math.Sqrt(Math.Max(0.0, disc));
        double sigmaMaxSq = (f + sqrtDisc) / 2.0;
        double sigmaMax = Math.Sqrt(sigmaMaxSq);
        return sigmaMax == 0.0 ? 0.0 : detAbs / sigmaMax;
    }

    /// <summary>Build the Householder reflector H = I − β·v·vᴴ that maps the sub-column x[start..end) to
    /// α·e_start (|α| = ‖x[start..end)‖), applies it to x in place (x[start] ← α, x[i&gt;start] ← 0), and
    /// returns α together with the reflector (v, β) for reuse on the other column. The phase choice
    /// α = −(x0/|x0|)·‖x‖ maximizes |v[start]|, avoiding cancellation.</summary>
    private static Complex HouseholderColumn(Complex[] x, int start, int end, out Complex[]? v, out double beta)
    {
        double xnorm = 0;
        for (int i = start; i < end; i++) xnorm += x[i].Real * x[i].Real + x[i].Imaginary * x[i].Imaginary;
        xnorm = Math.Sqrt(xnorm);
        if (xnorm == 0.0) { v = null; beta = 0.0; return Complex.Zero; }

        Complex x0 = x[start];
        Complex phase = x0.Magnitude == 0.0 ? Complex.One : x0 / x0.Magnitude;
        Complex alpha = -phase * xnorm;

        v = new Complex[end];
        v[start] = x0 - alpha;
        for (int i = start + 1; i < end; i++) v[i] = x[i];
        double vnorm2 = 0;
        for (int i = start; i < end; i++) vnorm2 += v[i].Real * v[i].Real + v[i].Imaginary * v[i].Imaginary;
        beta = vnorm2 == 0.0 ? 0.0 : 2.0 / vnorm2;

        // apply to x itself: x <- α e_start
        x[start] = alpha;
        for (int i = start + 1; i < end; i++) x[i] = Complex.Zero;
        return alpha;
    }

    /// <summary>Apply the reflector H = I − β·v·vᴴ to y[start..end) in place: y ← y − β·(vᴴy)·v.</summary>
    private static void ApplyReflector(Complex[]? v, double beta, Complex[] y, int start, int end)
    {
        if (v is null || beta == 0.0) return;
        Complex s = Complex.Zero;
        for (int i = start; i < end; i++) s += Complex.Conjugate(v[i]) * y[i];
        Complex factor = beta * s;
        for (int i = start; i < end; i++) y[i] -= factor * v[i];
    }

    // -----------------------------------------------------------------------------------------------
    // the (1,2) witness plane + the full member pipeline
    // -----------------------------------------------------------------------------------------------

    /// <summary>The near-defective 2-plane at the seed's (1,2) block, expanded to the FULL Configs basis:
    /// refine the seed ONCE (<see cref="SectorShellCensus.RefineSeed"/>), build the (1,2) sector in the
    /// seed's R-parity, subtract λ_A, and take the two smallest right singular vectors of the SHIFTED
    /// sector [N3] (the invariant 2-subspace is well-conditioned though the individual eigenvectors are
    /// near-parallel), then expand each to the full block. Returns the two full plane vectors together
    /// with the refined λ_A, q, and pair gap.</summary>
    public static (Complex[] V1, Complex[] V2, Complex Lambda, Complex Q, double PairGap) WitnessPlane(RealSeed seed)
    {
        var (qRefined, lambda, pairGap) = SectorShellCensus.RefineSeed(seed);
        var q = new Complex(qRefined, 0);
        bool odd = seed.RParity < 0;

        var (aSec, dSec) = WeightCoherenceBlock.BuildReflectionSectorColumnMajor(seed.N, 1, 2, q, odd);
        var m = Matrix<Complex>.Build.Dense(dSec, dSec);
        for (int c = 0; c < dSec; c++)
            for (int r = 0; r < dSec; r++)
                m[r, c] = aSec[(long)c * dSec + r] - (r == c ? lambda : Complex.Zero);

        var svd = m.Svd(computeVectors: true);
        var vt = svd.VT;                                                    // VT = Vᴴ; right sing. vec i = conj(row i)
        var v1sec = new Complex[dSec];
        var v2sec = new Complex[dSec];
        for (int i = 0; i < dSec; i++)
        {
            v1sec[i] = Complex.Conjugate(vt[dSec - 1, i]);                  // smallest singular value
            v2sec[i] = Complex.Conjugate(vt[dSec - 2, i]);                  // second smallest
        }

        var v1 = ExpandSectorToFull(seed.N, 1, 2, odd, v1sec);
        var v2 = ExpandSectorToFull(seed.N, 1, 2, odd, v2sec);
        return (v1, v2, lambda, q, pairGap);
    }

    /// <summary>Expand an R-parity sector vector to the full (wKet,wBra) Configs basis, mirroring
    /// <see cref="WeightCoherenceBlock.BuildReflectionSectorColumnMajor"/>'s basis: a reflection fixed
    /// point carries weight 1, a 2-cycle representative t &lt; Rt maps to (e_t ± e_{Rt})/√2 (sign +
    /// even, − odd). The reps ordering (fixed points first for the even sector, then 2-cycle reps in
    /// increasing index) is transcribed verbatim so the SVD sector components line up.</summary>
    public static Complex[] ExpandSectorToFull(int n, int wKet, int wBra, bool odd, Complex[] sectorVec)
    {
        if (sectorVec is null) throw new ArgumentNullException(nameof(sectorVec));
        var perm = WeightCoherenceBlock.ReflectionPermutation(n, wKet, wBra);
        var reps = new List<int>();
        for (int i = 0; i < perm.Length; i++)
        {
            if (perm[i] == i) { if (!odd) reps.Add(i); }
            else if (i < perm[i]) reps.Add(i);
        }
        if (sectorVec.Length != reps.Count)
            throw new ArgumentException($"sector vector length {sectorVec.Length} != sector dim {reps.Count}", nameof(sectorVec));

        double sSign = odd ? -1.0 : 1.0;
        double inv2 = 1.0 / Math.Sqrt(2.0);
        var full = new Complex[perm.Length];
        for (int r = 0; r < reps.Count; r++)
        {
            int t = reps[r];
            if (perm[t] == t) full[t] += sectorVec[r];
            else { full[t] += sectorVec[r] * inv2; full[perm[t]] += sectorVec[r] * (sSign * inv2); }
        }
        return full;
    }

    /// <summary>The full member pipeline for a seed: refine once [N9], build the (1,2) witness plane in
    /// the seed's R-parity [N3], climb the |p−q̃|=1 band with <see cref="ApplyW"/> (re-orthonormalizing
    /// the transported pair after EVERY step [N3]), and at each band member (p,p+1) in
    /// <see cref="SectorShellCensus.ExpectedMembers"/> read a from-above σ_min bound at λ_A and, via
    /// <see cref="ApplyFoldGauge"/>, at the folded μ = −λ_A−2N member. Keys match ExpectedMembers
    /// exactly; parity bookkeeping follows the alternation law [F7]. The value at each key is the
    /// tightest (minimum) bound found for that cell.</summary>
    public static IReadOnlyDictionary<(int P, int W, string Shift), (double Bound, bool CarriedOdd)>
        MemberUpperBounds(RealSeed seed, Action<string>? log = null)
    {
        int n = seed.N;
        var (v1, v2, lambda, q, _) = WitnessPlane(seed);
        Complex mu = -lambda - 2.0 * n;
        bool seedOdd = seed.RParity < 0;

        var result = new Dictionary<(int P, int W, string Shift), (double Bound, bool CarriedOdd)>();
        Complex[] c1 = v1, c2 = v2;                                         // the plane, currently at (1,2)

        // Climb only to the last IN-FUNDAMENTAL-DOMAIN band block (p ≤ (N−1)/2, i.e. 2p+1 ≤ N): every
        // recorded band member sits there (band ∈ FD ⟺ 2p+1 ≤ N), and every μ rep(p,N−p−1) has its
        // smaller index ≤ (N−1)/2, so its band source is already in this range. This never climbs W
        // PAST half filling, where W is rank-deficient (§7: e.g. (3,4)→(4,5) shrinks the dimension);
        // the upper rungs are reached by the fold, exactly as the corollary prescribes.
        int pMax = (n - 1) / 2;
        for (int p = 1; p <= pMax; p++)
        {
            int w = p + 1;
            bool carriedOdd = seedOdd ^ (((p - 1) & 1) == 1);              // alternation law [F7]

            // band member (p,p+1) at λ_A (always in the fundamental-domain strip for p ≤ pMax)
            if (BlockLattice.InFundamentalDomain(n, p, w))
            {
                var csr = WeightCoherenceSectorCsr.BuildFull(n, p, w, q);
                double bound = TwoPlaneResidual(csr, lambda, c1, c2);
                Merge(result, (p, w, "lambdaA"), bound, carriedOdd);
                log?.Invoke($"band ({p},{w}) x lambdaA: bound={bound:E3} carriedOdd={carriedOdd}");
            }

            // fold to the μ member (p, N-w); the fold preserves the carried parity
            {
                int fw = n - w;                                            // = n - p - 1
                var f1 = ApplyFoldGauge(n, p, w, c1);
                var f2 = ApplyFoldGauge(n, p, w, c2);
                var csrFold = WeightCoherenceSectorCsr.BuildFull(n, p, fw, q);
                double boundMu = TwoPlaneResidual(csrFold, mu, f1, f2);
                var rep = p <= fw ? (p, fw) : (fw, p);                     // transpose-normalized FD rep
                Merge(result, (rep.Item1, rep.Item2, "mu"), boundMu, carriedOdd);
                log?.Invoke($"fold ({p},{w})->({p},{fw}) x mu: bound={boundMu:E3} rep=({rep.Item1},{rep.Item2}) carriedOdd={carriedOdd}");
            }

            // climb to (p+1,p+2), re-orthonormalizing the transported pair [N3]
            if (p < pMax)
            {
                c1 = ApplyW(n, p, w, c1);
                c2 = ApplyW(n, p, w, c2);
                Orthonormalize(c1, c2);
            }
        }
        return result;
    }

    private static void Merge(Dictionary<(int, int, string), (double Bound, bool CarriedOdd)> d,
        (int, int, string) key, double bound, bool carriedOdd)
    {
        if (d.TryGetValue(key, out var cur))
        { if (bound < cur.Bound) d[key] = (bound, carriedOdd); }           // keep the tightest upper bound
        else d[key] = (bound, carriedOdd);
    }

    // -----------------------------------------------------------------------------------------------
    // small complex-vector helpers
    // -----------------------------------------------------------------------------------------------

    private static Dictionary<int, int> Index(List<int> configs)
    {
        var index = new Dictionary<int, int>(configs.Count);
        for (int i = 0; i < configs.Count; i++) index[configs[i]] = i;
        return index;
    }

    /// <summary>In-place 2-column modified Gram-Schmidt orthonormalization of the transported pair.</summary>
    private static void Orthonormalize(Complex[] a, Complex[] b)
    {
        double na = Norm(a);
        if (na > 0) Scale(a, 1.0 / na);
        Complex proj = Inner(a, b);
        for (int i = 0; i < b.Length; i++) b[i] -= proj * a[i];
        double nb = Norm(b);
        if (nb > 0) Scale(b, 1.0 / nb);
    }

    private static Complex Inner(Complex[] a, Complex[] b)                  // <a,b> = Σ conj(a)·b
    {
        Complex s = Complex.Zero;
        for (int i = 0; i < a.Length; i++) s += Complex.Conjugate(a[i]) * b[i];
        return s;
    }

    private static double Norm(Complex[] v)
    {
        double s = 0;
        foreach (var z in v) s += z.Real * z.Real + z.Imaginary * z.Imaginary;
        return Math.Sqrt(s);
    }

    private static void Scale(Complex[] v, double f)
    {
        for (int i = 0; i < v.Length; i++) v[i] *= f;
    }
}
