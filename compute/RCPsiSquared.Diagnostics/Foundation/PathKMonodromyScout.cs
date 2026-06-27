using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Generalises the path-3 octic monodromy (Gal(F_8)=S_8, --root galoismonodromy) off path-3 to
/// arbitrary path-k. The (SE,DE) S₂-symmetric block 2M(q)=A+q·C is built for any nBlock=k+1 via
/// <see cref="F89PathKSeDeBlock"/>; the AT-locked half is identified exactly by matching the block spectrum
/// at the base q0 against the roots of the reconstructed AT factor (<see cref="F89AtFactorReconstruction"/>,
/// the rate-confined invariant subspace, since the naive Slater rule fails from path-5); the remaining
/// residual strands (degree F_d = 8/18/32/53 for k=3/4/5/6) carry the H_B-mixed Galois group. As q loops the
/// complex plane the residual strands braid; lassoing every branch point (EP) from a common base and
/// union-finding the transposition graph reconstructs the monodromy = Galois group from below (S_d ⟺ the
/// residual graph is connected). A SCOUT: it validates against the known Gal = S_8 (path-3) / S_18 (path-4)
/// algebraic certificates, and maps the entirely-uncomputed path-4+ branch geometry (EP count, locations,
/// any diabolic point). The within-block σ_T fold (zeros/twins) is deliberately NOT carried over: it is
/// N=4-only (foldcross), so the cross-block "± is the road" layer is a separate, later step.
///
/// <para>Reuses the generic <see cref="Monodromy"/> tracker; the path-3 witness <see cref="GaloisMonodromyWitness"/>
/// stays the validated N=4 special case.</para></summary>
public static class PathKMonodromyScout
{
    private static Complex[,] ToComplex(RCPsiSquared.Core.Numerics.GaussianInteger[,] g)
    {
        int d = g.GetLength(0);
        var c = new Complex[d, d];
        for (int r = 0; r < d; r++)
            for (int s = 0; s < d; s++)
                c[r, s] = new Complex((double)g[r, s].Re, (double)g[r, s].Im);
        return c;
    }

    // 2M(q) = A + q·C, extracted once at q0=0 (=A) and q0=1 (C = 2M(1)−A). Eigenvalues of (A+qC)/2 are λ_k.
    public static (Complex[,] A, Complex[,] C) BuildLinear(int nBlock)
    {
        var a = ToComplex(F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0: 0, nBlock: nBlock));
        var m1 = ToComplex(F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0: 1, nBlock: nBlock));
        int d = a.GetLength(0);
        var c = new Complex[d, d];
        for (int r = 0; r < d; r++)
            for (int s = 0; s < d; s++)
                c[r, s] = m1[r, s] - a[r, s];
        return (a, c);
    }

    // all block roots λ_k at complex q.
    public static Complex[] AllRootsAt(Complex[,] a, Complex[,] c, Complex q)
    {
        int d = a.GetLength(0);
        var m = Matrix<Complex>.Build.Dense(d, d, (r, s) => a[r, s] + q * c[r, s]);
        return m.Evd().EigenValues.Select(v => v / 2).ToArray();
    }

    // the AT-locked eigenvalues λ_AT at integer q0, from the exact reconstructed AT factor (roots 2·λ_AT).
    private static Complex[] AtRootsAt(int k, int q0)
    {
        var coeffs = F89AtFactorReconstruction.ForPathK(k);             // monic, lowest-first, ×2 scale
        // the AT factor is built at q0=2 inside ForPathK; rescale roots by the linear-in-q structure is not
        // needed — we evaluate identification at the SAME q0=2 base the factor is built for.
        int d = coeffs.Length - 1;
        var comp = Matrix<Complex>.Build.Dense(d, d);                   // companion matrix of the monic poly
        for (int i = 0; i < d; i++)
            comp[i, d - 1] = -new Complex((double)coeffs[i].Re, (double)coeffs[i].Im);
        for (int i = 1; i < d; i++) comp[i, i - 1] = Complex.One;
        return comp.Evd().EigenValues.Select(v => v / 2).ToArray();    // 2·λ_AT → λ_AT
    }

    /// <summary>The residual (H_B-mixed, non-AT) start-indices among the block roots at q0, identified by
    /// removing the strands nearest the exact AT eigenvalues. Returns (residualIndices, atIndices).</summary>
    public static (int[] residual, int[] at) ResidualIndices(int k, Complex[] r0)
    {
        var atRoots = AtRootsAt(k, 2);
        var atIdx = new HashSet<int>();
        foreach (var t in atRoots)
        {
            int best = -1; double bd = double.PositiveInfinity;
            for (int i = 0; i < r0.Length; i++)
                if (!atIdx.Contains(i)) { double dd = (r0[i] - t).Magnitude; if (dd < bd) { bd = dd; best = i; } }
            atIdx.Add(best);
        }
        var residual = Enumerable.Range(0, r0.Length).Where(i => !atIdx.Contains(i)).ToArray();
        return (residual, atIdx.OrderBy(x => x).ToArray());
    }

    /// <summary>The F_d residual roots λ_k(q) at q for path-k: <see cref="AllRootsAt"/> with the AT-locked
    /// strands removed. The path-k analogue of <see cref="GaloisMonodromyWitness.OcticRootsAt"/>; reproduces
    /// the path-3 octic exactly at q=2 (the trusted cross-check).
    /// <para><b>SCOPE (R-7): valid ONLY at q≈2.</b> <see cref="ResidualIndices"/> matches against
    /// <c>AtRootsAt(k, 2)</c> — the AT factor reconstructed at a FIXED q0=2 (the naive Slater rule that would
    /// give analytic q-scaling fails from path-5). The block's true AT eigenvalues move with q (Im ∝ J=qγ),
    /// so away from q=2 the nearest-match mis-labels which strands are AT. Use this for the single-q
    /// cross-check only; the diabolic gap-scan tracks the residual SET from q=2 by continuity instead (the
    /// residual set is monodromy-invariant: F_AT·F_residual are distinct factors). See FindDiabolics.</para></summary>
    public static Complex[] ResidualRootsAt(int k, Complex q)
    {
        var (a, c) = BuildLinear(k + 1);                 // nBlock = k+1
        var all = AllRootsAt(a, c, q);
        var (residual, _) = ResidualIndices(k, all);     // tuple (residual, at); AT-locked stripped via ForPathK
        return residual.Select(i => all[i]).ToArray();
    }

    /// <summary>The σ_T classification of the residual strands under the GLOBAL palindrome fold
    /// λ ↦ −λ̄ − 2σ (σ = nBlock = −N centre). For each residual strand i, FoldPartner[i] = j if its
    /// fold-image is residual strand j (within-block): i==j ⟹ an on-fold "zero" (self-mirror, Re λ ≈ −σ),
    /// i≠j ⟹ a within-block "twin" (a ± mode). FoldPartner[i] = −1 ⟹ CROSS-block: the fold-image is not in
    /// the (SE,DE) residual at all, so the strand's mirror partner lives in the (SE,w_{N−2}) block (foldcross).
    /// At path-3 (partner=self, N=4) this reproduces the within-block σ_T: 4 zeros + 2 twin-pairs. At path-k≥4
    /// every strand is cross-block (0 zeros), since the (SE,DE) block no longer self-folds.</summary>
    public static int[] FoldPartners(Complex[] residualRoots, int nBlock, double tol = 1e-4)
    {
        double sigma = nBlock;
        var part = new int[residualRoots.Length];
        for (int i = 0; i < residualRoots.Length; i++)
        {
            var img = new Complex(-residualRoots[i].Real - 2 * sigma, residualRoots[i].Imaginary);
            int best = -1; double bd = double.PositiveInfinity;
            for (int j = 0; j < residualRoots.Length; j++)
            {
                double d = (residualRoots[j] - img).Magnitude;
                if (d < bd) { bd = d; best = j; }
            }
            part[i] = bd < tol ? best : -1;            // −1 = cross-block (fold-image outside the residual)
        }
        return part;
    }

    private static int Moved(int[] perm) => Enumerable.Range(0, perm.Length).Count(i => perm[i] != i);

    // the cycles (length ≥ 2) of perm restricted to the residual strands, as lists of residual positions
    // (0..rd−1). Called only when every moved index is residual, so each cycle stays within the residual set.
    private static List<int[]> ResidualCycles(int[] perm, Dictionary<int, int> pos)
    {
        var cycles = new List<int[]>();
        var seen = new HashSet<int>();
        foreach (var idx in pos.Keys)
        {
            if (seen.Contains(idx) || perm[idx] == idx) continue;
            var cyc = new List<int>();
            int cur = idx;
            while (!seen.Contains(cur)) { seen.Add(cur); cyc.Add(pos[cur]); cur = perm[cur]; }
            if (cyc.Count >= 2) cycles.Add(cyc.ToArray());
        }
        return cycles;
    }

    // detour lasso bridging over the q=0 super-branch (the q^(big) factor), as in the path-3 witness.
    private static Complex[] DetourLasso(Complex q0, Complex ep, double radius)
    {
        double sgn = ep.Imaginary >= 0 ? 1.0 : -1.0;
        double h = 1.8 * sgn;
        Complex up0 = new(q0.Real, h), upE = new(ep.Real, h);
        Complex enter = ep + new Complex(0, radius * sgn);
        var pts = new List<Complex> { q0 };
        void Line(Complex a, Complex b)
        {
            int n = Math.Max(40, (int)(120 * (b - a).Magnitude));
            for (int kk = 1; kk <= n; kk++) pts.Add(a + (b - a) * ((double)kk / n));
        }
        Line(q0, up0); Line(up0, upE); Line(upE, enter);
        double th0 = Math.Atan2(enter.Imaginary - ep.Imaginary, enter.Real - ep.Real);
        const int nC = 240;
        for (int kk = 1; kk <= nC; kk++) { double th = th0 + 2 * Math.PI * kk / nC; pts.Add(ep + radius * new Complex(Math.Cos(th), Math.Sin(th))); }
        Line(enter, upE); Line(upE, up0); Line(up0, q0);
        return pts.ToArray();
    }

    // find branch points by the topological monodromy-cell test (an EP makes a small enclosing loop braid).
    public static List<Complex> FindBranchPoints(Func<Complex, Complex[]> roots,
        double reLo, double reHi, double imLo, double imHi, double cell, double mask = 0.20)
    {
        int nRe = Math.Max(1, (int)((reHi - reLo) / cell)), nIm = Math.Max(1, (int)((imHi - imLo) / cell));
        var fires = new bool[nRe, nIm];
        Parallel.For(0, nRe, ir =>
        {
            double re = reLo + ir * cell + cell / 2;
            for (int ii = 0; ii < nIm; ii++)
            {
                var c = new Complex(re, imLo + ii * cell + cell / 2);
                if (c.Magnitude >= mask && Moved(Monodromy.Permutation(roots, c, cell * 0.62, 80)) > 0)
                    fires[ir, ii] = true;
            }
        });

        var found = new List<Complex>();
        for (int ir = 0; ir < nRe; ir++)
            for (int ii = 0; ii < nIm; ii++)
            {
                if (!fires[ir, ii]) continue;
                var c0 = new Complex(reLo + ir * cell + cell / 2, imLo + ii * cell + cell / 2);
                if (found.Any(b => (b - c0).Magnitude < 1.5 * cell)) continue;
                var c = QuarterRefine(roots, c0, cell);
                if (found.Any(b => (b - c).Magnitude < 0.02)) continue;
                if (Moved(Monodromy.Permutation(roots, c, 0.6 * cell, 200)) > 0) found.Add(c);
            }
        return found;
    }

    private static Complex QuarterRefine(Func<Complex, Complex[]> roots, Complex c, double cell)
    {
        double half = cell / 2;
        for (int it = 0; it < 9 && half > 3e-4; it++)
        {
            double q = half / 2;
            Complex pick = c; bool found = false;
            foreach (var off in new[] { new Complex(q, q), new Complex(-q, q), new Complex(q, -q), new Complex(-q, -q) })
                if (Moved(Monodromy.Permutation(roots, c + off, half * 0.95, 100)) > 0) { pick = c + off; found = true; break; }
            c = pick; half = q;
            if (!found) { /* EP near centre, keep shrinking */ }
        }
        return c;
    }

    // ---- the diabolic hunt (the N=4→N=5 forward edge): gap-minima that DON'T braid (identity loop) ----

    /// <summary>A coalescence found by <see cref="FindDiabolics"/>. QValue = q* where two strands meet,
    /// MergeLambda = their midpoint λ_d; IsSemisimple = EpCharacter reads Diabolic (geo=alg, departure≈0,
    /// NOT a Jordan block); LoopIsIdentity = a small q-loop about q* braids nothing (necessary-not-sufficient
    /// for diabolic, R-3); Gap = refined min gap; GapScalingExponent = the gap's power-law in |q−q*|
    /// (≈1 linear ⟹ two sheets crossing ⟹ diabolic; ≈½ ⟹ defective √-branch); PairIsResidual = both
    /// coalescing strands are H_B-mixed residual (not AT-locked), tested by continuity from q=2.</summary>
    public sealed record DiabolicPoint(
        Complex QValue, Complex MergeLambda, bool IsSemisimple, bool LoopIsIdentity,
        double Gap, double GapScalingExponent, bool PairIsResidual);

    private static double MinGap(Complex[] r)
    {
        double m = double.PositiveInfinity;
        for (int i = 0; i < r.Length; i++)
            for (int j = i + 1; j < r.Length; j++)
                m = Math.Min(m, (r[i] - r[j]).Magnitude);
        return m;
    }

    private static Matrix<Complex> BlockAt(Complex[,] a, Complex[,] c, Complex q)
    {
        int d = a.GetLength(0);
        return Matrix<Complex>.Build.Dense(d, d, (r, s) => (a[r, s] + q * c[r, s]) / 2);
    }

    // gradient-free descent on the min-gap field toward a local coalescence q*.
    private static Complex GapRefine(Func<Complex, Complex[]> roots, Complex c, double cell)
    {
        double half = cell / 2, curGap = MinGap(roots(c));
        for (int it = 0; it < 40 && half > 1e-8; it++)
        {
            Complex best = c; double bestGap = curGap; bool moved = false;
            foreach (var off in new[]
            {
                new Complex(half, 0), new Complex(-half, 0), new Complex(0, half), new Complex(0, -half),
                new Complex(half, half), new Complex(half, -half), new Complex(-half, half), new Complex(-half, -half),
            })
            {
                double g = MinGap(roots(c + off));
                if (g < bestGap) { bestGap = g; best = c + off; moved = true; }
            }
            if (moved) { c = best; curGap = bestGap; } else half /= 2;
        }
        return c;
    }

    // classify the coalescing pair at q* by tracking the FULL block from q0=2 by continuity (R-4b/R-7): the
    // residual SET is monodromy-invariant (F_AT·F_residual are distinct factors, AT strands single-valued),
    // so SET membership survives any label permutations on the path. Returns whether the pair is
    // residual-residual, the merge λ (midpoint), and a λ-radius that encloses ONLY the pair.
    private static (bool pairIsResidual, Complex mergeLambda, double charRadius) ClassifyPair(
        Func<Complex, Complex[]> roots, Complex q0, Complex qd, HashSet<int> residualSet, int steps = 400)
    {
        var path = new Complex[steps + 1];
        for (int s = 0; s <= steps; s++) path[s] = q0 + (qd - q0) * ((double)s / steps);
        var traj = Monodromy.TrajectoryAlongPath(roots, path);
        var final = traj[steps];                       // final[startLabel] = position at q*

        int ai = 0, bi = 1; double best = double.PositiveInfinity;
        for (int i = 0; i < final.Length; i++)
            for (int j = i + 1; j < final.Length; j++)
            {
                double g = (final[i] - final[j]).Magnitude;
                if (g < best) { best = g; ai = i; bi = j; }
            }
        bool pairResidual = residualSet.Contains(ai) && residualSet.Contains(bi);
        var mid = (final[ai] + final[bi]) / 2;
        double distOther = double.PositiveInfinity;     // nearest non-pair strand → enclosing radius
        for (int i = 0; i < final.Length; i++)
            if (i != ai && i != bi) distOther = Math.Min(distOther, (final[i] - mid).Magnitude);
        double radius = Math.Clamp(0.4 * distOther, 0.01, 0.5);
        return (pairResidual, mid, radius);
    }

    // the gap's power-law exponent approaching q* along +Re q: gap ∝ |t|^p, p≈1 diabolic / p≈½ defective.
    private static double GapScalingExponent(Func<Complex, Complex[]> roots, Complex qd)
    {
        var ts = new[] { 0.02, 0.01, 0.005, 0.0025, 0.00125 };
        var xs = new List<double>(); var ys = new List<double>();
        foreach (var t in ts)
        {
            double g = MinGap(roots(qd + new Complex(t, 0)));
            if (g > 1e-13) { xs.Add(Math.Log(t)); ys.Add(Math.Log(g)); }
        }
        if (xs.Count < 2) return double.NaN;
        double mx = xs.Average(), my = ys.Average(), num = 0, den = 0;
        for (int i = 0; i < xs.Count; i++) { num += (xs[i] - mx) * (ys[i] - my); den += (xs[i] - mx) * (xs[i] - mx); }
        return den == 0 ? double.NaN : num / den;
    }

    /// <summary>Find the path-k residual's DIABOLIC points: the gap-minima that pkmono's EP finder discards
    /// (their monodromy loop is the identity). Scans the q-region for full-block min-gap dips, refines each,
    /// and classifies it from below — pair-is-residual (continuity from q=2), identity-loop, gap-scaling
    /// exponent, and the EpCharacter verdict (Diabolic vs Defective). The N=4→N=5 forward-edge instrument:
    /// at path-3 it must re-find q_EP at λ=−4+iJ·2; at path-k≥4 it answers whether a within-block diabolic
    /// survives once the self-fold is gone (R-1/R-2 predict it generically does NOT). The q=0 super-branch is
    /// masked per-k (default 0.20, NOT the octic's path-3 0.18). Returns ALL coalescences with refined gap
    /// below <paramref name="gapTol"/>; the caller filters IsSemisimple for genuine diabolics.</summary>
    public static List<DiabolicPoint> FindDiabolics(int k, double reLo, double reHi, double imLo, double imHi,
        double cell, double mask = 0.20, double gapTol = 1e-3, double loopRadius = 0.02)
    {
        var (a, c) = BuildLinear(k + 1);
        Func<Complex, Complex[]> roots = qq => AllRootsAt(a, c, qq);
        var q0 = new Complex(2, 0);
        var (residual, _) = ResidualIndices(k, roots(q0));
        var residualSet = new HashSet<int>(residual);

        int nRe = Math.Max(1, (int)((reHi - reLo) / cell)), nIm = Math.Max(1, (int)((imHi - imLo) / cell));
        var gap = new double[nRe, nIm];
        Parallel.For(0, nRe, ir =>
        {
            double re = reLo + ir * cell + cell / 2;
            for (int ii = 0; ii < nIm; ii++)
            {
                var q = new Complex(re, imLo + ii * cell + cell / 2);
                gap[ir, ii] = q.Magnitude < mask ? double.NaN : MinGap(roots(q));
            }
        });

        // candidate cells: a strict local minimum of the (smooth) gap field — a coalescence is a local min
        // that GapRefine drives toward 0. NO absolute threshold (a steep-slope diabolic has gap ~0.07 even at
        // the nearest cell, far above any cell-based cutoff; the q_EP nearest-cell gap is ~11.6·0.006); the
        // gapTol filter after refine rejects the avoided crossings. A smooth field has few local minima, so
        // this does not over-seed; the expensive classify runs only AFTER the gapTol filter.
        var seeds = new List<Complex>();
        for (int ir = 0; ir < nRe; ir++)
            for (int ii = 0; ii < nIm; ii++)
            {
                double g = gap[ir, ii];
                if (double.IsNaN(g)) continue;
                bool isMin = true;
                for (int dr = -1; dr <= 1 && isMin; dr++)
                    for (int di = -1; di <= 1; di++)
                    {
                        if (dr == 0 && di == 0) continue;
                        int jr = ir + dr, ji = ii + di;
                        if (jr < 0 || jr >= nRe || ji < 0 || ji >= nIm) continue;
                        if (!double.IsNaN(gap[jr, ji]) && gap[jr, ji] < g) { isMin = false; break; }
                    }
                if (isMin) seeds.Add(new Complex(reLo + ir * cell + cell / 2, imLo + ii * cell + cell / 2));
            }

        var result = new List<DiabolicPoint>();
        foreach (var seed in seeds)
        {
            var qd = GapRefine(roots, seed, cell);
            if (qd.Magnitude < mask) continue;                                       // refined INTO the q=0 super-branch (the q^big pile-up where the block goes diagonal, all rates collapse onto -2/-6) - not a real coalescence
            if (result.Any(d => (d.QValue - qd).Magnitude < 2 * cell)) continue;     // dedupe
            double g = MinGap(roots(qd));
            if (g > gapTol) continue;                                                // an avoided crossing, not a coalescence

            var (pairResidual, mergeLambda, charRadius) = ClassifyPair(roots, q0, qd, residualSet);
            bool loopId = Moved(Monodromy.Permutation(roots, qd, loopRadius, 240)) == 0;
            double expo = GapScalingExponent(roots, qd);
            var reading = EpCharacter.Characterize(BlockAt(a, c, qd), mergeLambda, charRadius);
            bool semisimple = reading.Kind == EpCharacter.EpKind.Diabolic;
            result.Add(new DiabolicPoint(qd, mergeLambda, semisimple, loopId, g, expo, pairResidual));
        }
        return result;
    }

    /// <summary>The full EpCharacter reading of the path-k residual block at (q, λ): the load-bearing
    /// discriminant (R-3) that separates a genuine DIABOLIC (semisimple, geo=alg, departure≈0) from a
    /// DEFECTIVE EP (Jordan block, geo&lt;alg) and from an avoided crossing (an identity loop alone cannot).
    /// Reuses the Riesz-projector classifier on the full block M=(A+qC)/2; <paramref name="radius"/> must
    /// enclose only the cluster at λ. The standalone entry point behind FindDiabolics' IsSemisimple.</summary>
    public static EpCharacter.Reading CharacterizeAt(int k, Complex q, Complex lambda, double radius = 0.1)
    {
        var (a, c) = BuildLinear(k + 1);
        return EpCharacter.Characterize(BlockAt(a, c, q), lambda, radius);
    }

    /// <summary>True iff the path-k residual block has a SEMISIMPLE (diabolic) repeated eigenvalue at (q, λ):
    /// EpCharacter reads Diabolic (eigenvectors stay independent). False for a defective EP or a simple
    /// eigenvalue. See <see cref="CharacterizeAt"/>.</summary>
    public static bool IsSemisimpleAt(int k, Complex q, Complex lambda, double radius = 0.1)
        => CharacterizeAt(k, q, lambda, radius).Kind == EpCharacter.EpKind.Diabolic;

    public sealed record ScanResult(
        int K, int NBlock, int BlockDim, int ResidualDim, int AtDim, Complex Q0,
        Complex[] ResidualRoots, IReadOnlyList<(Complex Q, int A, int B, int[] MovedResidual)> Eps,
        IReadOnlyList<(int A, int B)> Edges, int Components, int Largest, int[] StrandComponent,
        int[] FoldPartner);

    /// <summary>Scan the q-region for the path-k residual's branch points, lasso each from q0, and union-find
    /// the transposition graph on the residual strands. Connected ⟺ Gal(F_d) = S_d reconstructed from below.</summary>
    public static ScanResult Scan(int k, double reLo, double reHi, double imLo, double imHi, double cell, Complex q0)
    {
        int nBlock = k + 1;
        var (a, c) = BuildLinear(nBlock);
        int blockDim = a.GetLength(0);
        Func<Complex, Complex[]> roots = qq => AllRootsAt(a, c, qq);

        var r0 = roots(q0);
        var (residual, at) = ResidualIndices(k, r0);
        int rd = residual.Length;
        var pos = new Dictionary<int, int>();
        for (int p = 0; p < rd; p++) pos[residual[p]] = p;             // residual start-index → 0..rd−1
        var residualRoots = residual.Select(i => r0[i]).ToArray();

        var eps = FindBranchPoints(roots, reLo, reHi, imLo, imHi, cell);
        var perms = new int[eps.Count][];
        Parallel.For(0, eps.Count, j =>
            perms[j] = Monodromy.PermutationAlongPath(roots, DetourLasso(q0, eps[j], radius: 0.02)));

        var epInfo = new List<(Complex, int, int, int[])>();
        var edges = new List<(int, int)>();
        var uf = Enumerable.Range(0, rd).ToArray();
        int Find(int x) { while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; }

        for (int j = 0; j < eps.Count; j++)
        {
            var perm = perms[j];
            var moved = Enumerable.Range(0, perm.Length).Where(i => perm[i] != i).ToList();
            var movedResidual = moved.Where(pos.ContainsKey).Select(i => pos[i]).OrderBy(x => x).ToArray();
            // Union within each cycle of perm restricted to the residual strands. A clean lasso nets one
            // transposition; a lasso that encloses a near-degenerate EP cluster nets a product of cycles, and
            // each cycle still certifies its strands share a monodromy orbit (so union within it). Only an
            // AT-strand graze (a tracking near-degeneracy, AT roots are single-valued and must not move) is
            // rejected. Denser path-k spectra hit such grazes more, so this is more robust than requiring a
            // single clean cycle, and it stays valid (re-validated: path-3 still reconstructs a clean S_8).
            if (moved.Count >= 2 && moved.All(pos.ContainsKey))
            {
                foreach (var cyc in ResidualCycles(perm, pos))
                    for (int m = 1; m < cyc.Length; m++)
                    {
                        edges.Add((cyc[0], cyc[m]));
                        uf[Find(cyc[0])] = Find(cyc[m]);
                    }
                epInfo.Add((eps[j], movedResidual.Length > 0 ? movedResidual[0] : -1,
                            movedResidual.Length > 0 ? movedResidual[^1] : -1, movedResidual));
            }
            else epInfo.Add((eps[j], -1, -1, movedResidual));
        }

        int comps = Enumerable.Range(0, rd).Select(Find).Distinct().Count();
        var sizes = new Dictionary<int, int>();
        for (int i = 0; i < rd; i++) { int rt = Find(i); sizes[rt] = sizes.GetValueOrDefault(rt) + 1; }
        return new ScanResult(k, nBlock, blockDim, rd, at.Length, q0, residualRoots, epInfo, edges,
            comps, sizes.Count == 0 ? 0 : sizes.Values.Max(), Enumerable.Range(0, rd).Select(Find).ToArray(),
            FoldPartners(residualRoots, nBlock));
    }
}
