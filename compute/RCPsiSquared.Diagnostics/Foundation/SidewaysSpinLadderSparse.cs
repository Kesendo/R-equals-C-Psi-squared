using System.Numerics;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The chain walk that carries the sideways ladder past the dense-EIGENSOLVER TIME wall to N=9
/// (the falsifiable next case of <c>SidewaysSpinLadderClaim</c>: ℓ = 3, predicted S⁺ transport norms
/// √6, √10, √12, √12, √10, √6, confirmed 2026-08-09). At N=9 the middle blocks are 10584² and 15876², too
/// big for dense zgeev in reasonable time, but the walk never needed full spectra: per source block it
/// needs ONE eigenpair near the known fold value −λ_A − 2N, which one dense LU of (L − σ) plus exact
/// inverse iteration delivers (<see cref="EigNearShift"/>). The wall crossed is TIME, not memory: the LU
/// still materialises one dense block (4.03 GB at 15876², a >2 GB array, so the host process needs
/// AllowVeryLargeObjects; the test project sets it, the Cli does not, which is why the live witness stays
/// at N ≤ 7). Only the CONTROL half and the ladder applies are genuinely sparse
/// (≤ 2(N−1)+1 nonzeros per L row, ≤ N per ladder column).
///
/// <para>What is checked and what is imported. CONTROL (both ladders intertwine exactly, per column,
/// compared to 0.0) is recomputed here sparsely, column by column: (L_t·M − M·L_s)e_j from ladder applies
/// and per-column sparse accumulation, the same case-1 comparison as the dense walk. The rung norms come
/// from dense zgeev below <see cref="DenseCutoff"/> and from LU shift-invert above it. SHAPE and COUNT
/// are NOT re-derived at N=9: the diamond membership at N=9 is census-verified in full (the step-3 σ_min
/// shell census, see <c>SpectatorIntertwinerClaim</c>), and the walk verifies per source that the fold
/// eigenvalue is present, so what this class adds at N=9 is the transport-norm measurement, not the
/// sector census.</para>
///
/// <para>Conventions are <see cref="SidewaysSpinLadderChain"/>'s exactly (PopcountStates ascending, flat
/// = pIdx·Mq + qIdx, γ = 1, hop amplitude J = 2q*, no JW sign in the hop, the sign lives in the ladders);
/// the CSR/column builders and the vector applies are gated entry-identical against the dense builders at
/// N=5, and the whole walk is diffed rung by rung against the dense walk at N=7 with the LU path forced,
/// before N=9 is read (<c>SidewaysSpinLadderSparseTests</c>).</para></summary>
public static class SidewaysSpinLadderSparse
{
    /// <summary>Blocks at or above this dimension take the LU shift-invert path (N=9: 3024-dim blocks
    /// stay dense at ~a minute, 10584/15876 go LU); below it, dense zgeev.</summary>
    public const int DenseCutoff = 3100;

    /// <summary>One rung of the sparse walk: the fold eigenpair of the source (λ within the derived
    /// tolerance of −λ_A − 2N), the S⁺ transport norm, the image's eigenvector residual on the target
    /// (relative), and which path produced the eigenvector.</summary>
    public readonly record struct SparseRung(
        int P, int Q, Complex Lambda, double Norm, double EigResidual,
        int SourceDim, int TargetDim, bool ViaInverseIteration);

    /// <summary>One chain of the sparse walk; Terminal/TerminalRatio/Survivors as in the dense walk
    /// (the terminal source block is small at every N walked, so the survivors control stays dense).</summary>
    public readonly record struct SparseChain(
        int Total, IReadOnlyList<SparseRung> Rungs, double Terminal, double TerminalRatio,
        int Survivors, int TerminalSourceDim, int TerminalTargetDim);

    public sealed record SparseWalkResult(
        int N, double QStar, double LambdaA, double Fold,
        double WorstControlResidual, double FoldTolerance,
        IReadOnlyList<double> CgPredicted, IReadOnlyList<SparseChain> Chains);

    /// <summary>The recorded seed per N: 5 and 7 delegate to <see cref="SidewaysSpinLadderChain.Seed"/>
    /// (the dense gate's loci); 9 is the first R-parity +1 census entry of
    /// <see cref="RealDefectiveSeeds"/>, (q*, λ_A) = (0.591760, −5.1206).</summary>
    public static (double QStar, double LambdaA) Seed(int n) => n == 9
        ? (0.591760, -5.1206)
        : SidewaysSpinLadderChain.Seed(n);

    // ------------------------------------------------------------------- sparse block construction

    private static int Binom(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        long r = 1;
        for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
        return (int)r;
    }

    private static bool Interior(int n, int p, int q) => p >= 1 && p <= n - 1 && q >= 1 && q <= n - 1;

    private static double Abs2(Complex x) => x.Real * x.Real + x.Imaginary * x.Imaginary;

    /// <summary>All nearest-neighbour hop partners of computational state a (big-endian sites, bit
    /// n−1−l = site l): a with the occupations of sites l, l+1 swapped, for each l where exactly one of
    /// the two is occupied. Amplitude is the caller's J; no JW sign (adjacent-site hop).</summary>
    private static IEnumerable<long> HopPartners(int n, long a)
    {
        for (int l = 0; l < n - 1; l++)
        {
            long m1 = 1L << (n - 1 - l);
            long m2 = 1L << (n - 2 - l);
            bool o1 = (a & m1) != 0, o2 = (a & m2) != 0;
            if (o1 != o2) yield return a ^ m1 ^ m2;
        }
    }

    /// <summary>The (p,q̃) block of L as CSR (row-major, columns ascending): diagonal −2·popcount(a⊕b),
    /// off-diagonals −iJ from the ket hop (a′,b) ← (a,b) and +iJ from the bra hop (a,b′) ← (a,b), read
    /// row-wise via H's symmetry (the entry in row (a,b) from column (a′,b) is −iJ for every hop
    /// neighbour a′ of a, likewise +iJ for b). Entry-identical to
    /// <see cref="SidewaysSpinLadderChain.BuildBlock"/>, gated at N=5.</summary>
    public static WeightCoherenceSectorCsr.Csr BuildCsr(int n, int p, int qt, double j)
    {
        var statesP = BlockBasis.PopcountStates(n, p);
        var statesQ = BlockBasis.PopcountStates(n, qt);
        var pIdx = new Dictionary<long, int>();
        for (int i = 0; i < statesP.Count; i++) pIdx[statesP[i]] = i;
        var qIdx = new Dictionary<long, int>();
        for (int i = 0; i < statesQ.Count; i++) qIdx[statesQ[i]] = i;
        int mq = statesQ.Count, dim = statesP.Count * mq;

        var rowPtr = new int[dim + 1];
        var cols = new List<int>(dim * (2 * n - 1));
        var vals = new List<Complex>(dim * (2 * n - 1));
        var rowCols = new List<(int Col, Complex Val)>(2 * n - 1);
        for (int i = 0; i < statesP.Count; i++)
            for (int k = 0; k < mq; k++)
            {
                long a = statesP[i], b = statesQ[k];
                rowCols.Clear();
                rowCols.Add((i * mq + k,
                    new Complex(-2.0 * System.Numerics.BitOperations.PopCount((ulong)(a ^ b)), 0)));
                foreach (long a2 in HopPartners(n, a))
                    rowCols.Add((pIdx[a2] * mq + k, new Complex(0, -j)));
                foreach (long b2 in HopPartners(n, b))
                    rowCols.Add((i * mq + qIdx[b2], new Complex(0, j)));
                rowCols.Sort((x, y) => x.Col.CompareTo(y.Col));
                foreach (var (c, v) in rowCols) { cols.Add(c); vals.Add(v); }
                rowPtr[i * mq + k + 1] = cols.Count;
            }
        return new WeightCoherenceSectorCsr.Csr(dim, rowPtr, cols.ToArray(), vals.ToArray());
    }

    /// <summary>Apply a ladder to a vector: y = M·x with M = S⁺ (phi: false) or Φ (phi: true), block
    /// (p,q̃) → (p+1, q̃±1), the same coefficients as <see cref="SpectatorIntertwiner.BuildSPlus"/> /
    /// <see cref="SpectatorIntertwiner.BuildW"/> (gated identical at N=5). y must have the target
    /// block's dimension; it is cleared first.</summary>
    public static void ApplyLadder(int n, int p, int qt, bool phi, Complex[] x, Complex[] y)
    {
        var inP = BlockBasis.PopcountStates(n, p);
        var inQ = BlockBasis.PopcountStates(n, qt);
        var outP = BlockBasis.PopcountStates(n, p + 1);
        var outQ = BlockBasis.PopcountStates(n, qt + (phi ? 1 : -1));
        var outPIdx = new Dictionary<long, int>();
        for (int i = 0; i < outP.Count; i++) outPIdx[outP[i]] = i;
        var outQIdx = new Dictionary<long, int>();
        for (int i = 0; i < outQ.Count; i++) outQIdx[outQ[i]] = i;
        Array.Clear(y);
        for (int i = 0; i < inP.Count; i++)
            for (int k = 0; k < inQ.Count; k++)
            {
                var xv = x[i * inQ.Count + k];
                if (xv == Complex.Zero) continue;
                long a = inP[i], b = inQ[k];
                for (int l = 0; l < n; l++)
                {
                    if (SpectatorIntertwiner.Occupied(n, a, l)) continue;
                    if (phi == SpectatorIntertwiner.Occupied(n, b, l)) continue;   // phi: need l∉b; S⁺: need l∈b
                    long nb = phi ? b | SpectatorIntertwiner.SiteMask(n, l) : b & ~SpectatorIntertwiner.SiteMask(n, l);
                    int sign = SpectatorIntertwiner.JwSign(n, a, l) * SpectatorIntertwiner.JwSign(n, b, l);
                    if (!phi && (l & 1) == 1) sign = -sign;
                    int row = outPIdx[a | SpectatorIntertwiner.SiteMask(n, l)] * outQ.Count + outQIdx[nb];
                    y[row] += sign * xv;
                }
            }
    }

    // -------------------------------------------------------------------------- CONTROL, sparsely

    /// <summary>The exact CONTROL residual of one rung, column by column: max-entry of
    /// L_t·(M·e_j) − M·(L_s·e_j) over all source columns j, every product sparse. The same case-1
    /// comparison to 0.0 as the dense walk (each entry's two sides draw the same multiset of addends);
    /// feasible at 15876² because e_j touches ≤ n ladder entries and ≤ 2(n−1)+1 block entries.</summary>
    public static double ControlResidual(int n, int p, int qt, bool phi, double j)
    {
        int sdim = Binom(n, p) * Binom(n, qt);
        int tdim = Binom(n, p + 1) * Binom(n, qt + (phi ? 1 : -1));
        if (sdim == 0 || tdim == 0) return 0.0;
        var lsCols = ToColumns(BuildCsr(n, p, qt, j));
        var ltCols = ToColumns(BuildCsr(n, p + 1, qt + (phi ? 1 : -1), j));

        var ej = new Complex[sdim];
        var mej = new Complex[tdim];
        var lsej = new Complex[sdim];
        var mlsej = new Complex[tdim];
        var acc = new Complex[tdim];
        double worst = 0.0;
        for (int col = 0; col < sdim; col++)
        {
            ej[col] = Complex.One;
            ApplyLadder(n, p, qt, phi, ej, mej);                     // M·e_j (≤ n entries)
            ej[col] = Complex.Zero;

            foreach (var (r, v) in lsCols[col]) lsej[r] = v;         // L_s·e_j = column col of L_s
            ApplyLadder(n, p, qt, phi, lsej, mlsej);                 // M·(L_s·e_j)
            foreach (var (r, _) in lsCols[col]) lsej[r] = Complex.Zero;

            // acc = L_t·(M·e_j) − M·(L_s·e_j), accumulated sparsely.
            for (int i = 0; i < tdim; i++) acc[i] = -mlsej[i];
            for (int i = 0; i < tdim; i++)
            {
                if (mej[i] == Complex.Zero) continue;
                foreach (var (r, v) in ltCols[i]) acc[r] += v * mej[i];
            }
            for (int i = 0; i < tdim; i++)
                if (acc[i] != Complex.Zero)
                    worst = Math.Max(worst, acc[i].Magnitude);
        }
        return worst;
    }

    private static List<(int Row, Complex Val)>[] ToColumns(WeightCoherenceSectorCsr.Csr m)
    {
        var cols = new List<(int, Complex)>[m.Dim];
        for (int i = 0; i < m.Dim; i++) cols[i] = new List<(int, Complex)>(4);
        for (int r = 0; r < m.Dim; r++)
            for (int k = m.RowPtr[r]; k < m.RowPtr[r + 1]; k++)
                cols[m.ColIdx[k]].Add((r, m.Values[k]));
        return cols;
    }

    // --------------------------------------------------------------- the fold eigenpair, two paths

    /// <summary>Shift-invert eigenpair via ONE dense LU of (L − σI) plus inverse iteration with EXACT
    /// solves (zgetrf + zgetrs through <see cref="MklDirect"/>), the same factorize-then-iterate pattern
    /// as <see cref="ShiftedSigmaMin"/>. Iterative solvers were tried first and are recorded as the
    /// wrong tool here: the identity-preconditioned BiCGStab breaks down or stalls at interior shifts of
    /// these non-normal blocks (measured at N=5 intermittently and at N=7 on every solve, both with σ on
    /// the fold and detuned), so the O(d³/3) factorization is bought once and every solve after it is
    /// exact. The shift stays FIXED at σ (partial pivoting handles the near-singularity; the solution
    /// amplifying the eigen-direction is the point), and the near-defective fold PAIR (gap ~5e-4 at the
    /// seeds), which mixes under any single-vector iteration, is separated afterwards by a 2×2
    /// Rayleigh-Ritz on the span of the last two iterates. Deterministic start (seeded). Returns the best
    /// eigenpair, its eigen residual ‖Lu − λu‖ (u unit), and the outer count.</summary>
    public static (Complex Lambda, Complex[] V, double EigResidual, int Outer) EigNearShift(
        WeightCoherenceSectorCsr.Csr block, Complex sigma, int seed,
        int maxOuter = 30, Action<string>? log = null)
    {
        int dim = block.Dim;
        var a = new Complex[(long)dim * dim];                        // column-major L − σI, LU'd in place
        for (int r = 0; r < dim; r++)
            for (int k = block.RowPtr[r]; k < block.RowPtr[r + 1]; k++)
                a[(long)block.ColIdx[k] * dim + r] = block.Values[k];
        for (int i = 0; i < dim; i++) a[(long)i * dim + i] -= sigma;
        var ipiv = new int[dim];
        MklDirect.LuFactorizeRaw(a, dim, ipiv);

        var rng = new Random(seed);
        var x = new Complex[dim];
        for (int i = 0; i < dim; i++) x[i] = new Complex(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
        Normalize(x);

        var lx = new Complex[dim];
        double best = double.PositiveInfinity;
        var bestV = (Complex[])x.Clone();
        Complex bestLambda = sigma;
        Complex[]? prevAccepted = null, lastAccepted = null;
        int sinceImprovement = 0;
        int outer;
        for (outer = 1; outer <= maxOuter; outer++)
        {
            MklDirect.LuSolveRaw(a, dim, ipiv, x, conjugateTranspose: false);   // x ← (L − σ)⁻¹ x, exact
            Normalize(x);
            prevAccepted = lastAccepted;
            lastAccepted = (Complex[])x.Clone();
            CsrOps.Multiply(block, x, lx);
            Complex num = Complex.Zero;
            for (int i = 0; i < dim; i++) num += Complex.Conjugate(x[i]) * lx[i];
            var lambda = num;                                        // x is unit
            double res = 0.0;
            for (int i = 0; i < dim; i++) res += Abs2(lx[i] - lambda * x[i]);
            res = Math.Sqrt(res);
            log?.Invoke($"outer {outer}: λ={lambda.Real:F6} eigres={res:E1}");
            if (res < best * 0.9)
            {
                best = res; Array.Copy(x, bestV, dim); bestLambda = lambda; sinceImprovement = 0;
            }
            else if (++sinceImprovement >= 4) break;                 // stalled at the achievable floor
            if (best < 1e-12) break;
        }

        // The pair split: the fold eigenvalue sits in a near-defective PAIR, the iterate converges into
        // the pair PLANE quickly but separates the members only at the slow rate |λ₁−σ|/|λ₂−σ| ≈ 1
        // (the measured λ oscillation between the members). Two successive iterates are two in-plane
        // mixes, so a 2×2 Rayleigh-Ritz on their span separates the members, leaving the out-of-plane
        // error. Kept only if it beats the single-vector iterate; measured at N=5: 1.3e-8 → 4.8e-10.
        if (lastAccepted is not null && prevAccepted is not null)
        {
            var b1 = lastAccepted;
            var b2 = (Complex[])prevAccepted.Clone();
            Complex ip12 = Complex.Zero;
            for (int i = 0; i < dim; i++) ip12 += Complex.Conjugate(b1[i]) * b2[i];
            for (int i = 0; i < dim; i++) b2[i] -= ip12 * b1[i];
            double n2 = 0.0;
            for (int i = 0; i < dim; i++) n2 += Abs2(b2[i]);
            if (n2 > 1e-24)
            {
                double inv = 1.0 / Math.Sqrt(n2);
                for (int i = 0; i < dim; i++) b2[i] *= inv;
                var lb1 = new Complex[dim];
                var lb2 = new Complex[dim];
                CsrOps.Multiply(block, b1, lb1);
                CsrOps.Multiply(block, b2, lb2);
                Complex h11 = Complex.Zero, h12 = Complex.Zero, h21 = Complex.Zero, h22 = Complex.Zero;
                for (int i = 0; i < dim; i++)
                {
                    h11 += Complex.Conjugate(b1[i]) * lb1[i];
                    h12 += Complex.Conjugate(b1[i]) * lb2[i];
                    h21 += Complex.Conjugate(b2[i]) * lb1[i];
                    h22 += Complex.Conjugate(b2[i]) * lb2[i];
                }
                Complex half = (h11 + h22) / 2.0;
                Complex disc = Complex.Sqrt(half * half - (h11 * h22 - h12 * h21));
                // Both Ritz values are computed; among those that BEAT the single-vector residual, the
                // one nearer σ wins. There is deliberately NO σ-distance acceptance guard: an
                // unconverged in-plane MIX has a Rayleigh value BETWEEN the members, i.e. often nearer
                // σ than either true member, so a guard of the form |θ−σ| ≤ |λ_best−σ| + const keeps
                // the mix and rejects the split (the N=9 first run did exactly that: it reported the
                // mix −12.87931 at residual 1.2e-6 where the members −12.878060 / −12.880829 sit at
                // 8.6e-14; caught by an independent review reproduction). Residual is the only judge.
                Complex candLam = bestLambda;
                Complex[]? candV = null;
                double candRes = best;
                foreach (var theta in new[] { half + disc, half - disc })
                {
                    Complex y1, y2;                                  // 2×2 eigenvector at θ
                    if ((h11 - theta).Magnitude > (h22 - theta).Magnitude) { y1 = h12; y2 = theta - h11; }
                    else { y1 = theta - h22; y2 = h21; }
                    double yn = Math.Sqrt(Abs2(y1) + Abs2(y2));
                    if (yn < 1e-300) continue;
                    y1 /= yn; y2 /= yn;
                    var u = new Complex[dim];
                    for (int i = 0; i < dim; i++) u[i] = y1 * b1[i] + y2 * b2[i];
                    var lu = new Complex[dim];
                    CsrOps.Multiply(block, u, lu);
                    Complex lam = Complex.Zero;
                    for (int i = 0; i < dim; i++) lam += Complex.Conjugate(u[i]) * lu[i];
                    double res2 = 0.0;
                    for (int i = 0; i < dim; i++) res2 += Abs2(lu[i] - lam * u[i]);
                    res2 = Math.Sqrt(res2);
                    bool beats = res2 < best;
                    bool nearer = candV is null || (lam - sigma).Magnitude < (candLam - sigma).Magnitude;
                    // prefer any true split member over the mix; among comparable members, the nearer one
                    if (beats && (candV is null || res2 < candRes * 1e-2 || (res2 < candRes * 1e2 && nearer)))
                    {
                        candLam = lam; candV = u; candRes = res2;
                    }
                }
                if (candV is not null)
                {
                    log?.Invoke($"pair split: θ={candLam.Real:F6} eigres={candRes:E1} (was {best:E1})");
                    best = candRes; bestV = candV; bestLambda = candLam;
                }
            }
        }
        return (bestLambda, bestV, best, outer);
    }

    private static void Normalize(Complex[] x)
    {
        double n2 = 0.0;
        for (int i = 0; i < x.Length; i++) n2 += Abs2(x[i]);
        double inv = 1.0 / Math.Sqrt(n2);
        for (int i = 0; i < x.Length; i++) x[i] *= inv;
    }

    /// <summary>Dense fold eigenpair (blocks below <see cref="DenseCutoff"/>): zgeev with vectors,
    /// argmin |λ − fold|, exactly the dense walk's selection.</summary>
    private static (Complex Lambda, Complex[] V) EigDense(WeightCoherenceSectorCsr.Csr block, Complex fold)
    {
        int dim = block.Dim;
        var a = new Complex[(long)dim * dim];
        for (int r = 0; r < dim; r++)
            for (int k = block.RowPtr[r]; k < block.RowPtr[r + 1]; k++)
                a[r + (long)block.ColIdx[k] * dim] = block.Values[k];
        var (w, vr) = MklDirect.EigenvaluesAndVectorsDirectRaw(a, dim);
        int jBest = 0;
        for (int j = 1; j < dim; j++)
            if ((w[j] - fold).Magnitude < (w[jBest] - fold).Magnitude) jBest = j;
        var v = new Complex[dim];
        for (int i = 0; i < dim; i++) v[i] = vr[jBest * dim + i];
        Normalize(v);
        return (w[jBest], v);
    }

    // ------------------------------------------------------------------------------ the sparse walk

    /// <summary>The sparse chain walk at the recorded seed for <paramref name="n"/>. Odd N with a
    /// recorded seed only (5, 7 for cross-validation against the dense walk; 9 is the target).
    /// <paramref name="denseCutoff"/> defaults to <see cref="DenseCutoff"/>; the N=7 gate lowers it to
    /// force the inverse-iteration path over blocks whose dense answer is known.</summary>
    public static SparseWalkResult Run(int n, Action<string>? log = null, int denseCutoff = DenseCutoff)
    {
        var (qStar, lamA) = Seed(n);
        double j = 2.0 * qStar;
        double fold = -lamA - 2.0 * n;
        int[] chains = { n - 1, n + 1 };

        // CONTROL, sparse and exact, both ladders on every rung of both chains.
        double worstControl = 0.0;
        foreach (bool phi in new[] { true, false })
        {
            int dq = phi ? 1 : -1;
            foreach (int total in chains)
                for (int p = 0; p <= n; p++)
                {
                    int q = total - p;
                    if (q < 0 || q > n || p + 1 > n || q + dq < 0 || q + dq > n) continue;
                    if (Binom(n, p) * Binom(n, q) == 0 || Binom(n, p + 1) * Binom(n, q + dq) == 0) continue;
                    double r = ControlResidual(n, p, q, phi, j);
                    worstControl = Math.Max(worstControl, r);
                    log?.Invoke($"control {(phi ? "phi" : "S+")} ({p},{q}): {r:E1}");
                }
        }

        // The walk: fold eigenpair per source, S⁺ transport, residual on the target.
        var cg = SidewaysSpinLadderChain.CgCoefficients(n);
        var chainResults = new List<SparseChain>();
        double bestOffset = double.PositiveInfinity;
        var pending = new List<(int Total, List<SparseRung> Rungs, double Terminal,
            int Survivors, int TermSrc, int TermTgt)>();
        foreach (int total in chains)
        {
            var walk = Enumerable.Range(0, n + 1).Where(p => total - p >= 0 && total - p <= n)
                                 .Select(p => (P: p, Q: total - p)).ToList();
            var rungs = new List<SparseRung>();
            double terminal = double.NaN;
            int survivors = 0, termSrc = 0, termTgt = 0;
            foreach (var (p, q) in walk.Take(walk.Count - 1))
            {
                int sdim = Binom(n, p) * Binom(n, q);
                int tdim = Binom(n, p + 1) * Binom(n, q - 1);
                if (sdim == 0 || tdim == 0) continue;
                if (!Interior(n, p, q)) continue;   // boundary sources never carry the fold (rate rung)
                var ls = BuildCsr(n, p, q, j);
                bool sparse = sdim >= denseCutoff;
                Complex lambda;
                Complex[] v;
                double eigRes;
                if (sparse)
                {
                    (lambda, v, eigRes, _) = EigNearShift(ls, fold, seed: 20260809 + 100 * p + q);
                    log?.Invoke($"({p},{q}) dim {sdim}: inverse iteration λ = {lambda.Real:F6}{lambda.Imaginary:+0.000000;-0.000000}i, eig residual {eigRes:E1}");
                }
                else
                {
                    (lambda, v) = EigDense(ls, fold);
                    log?.Invoke($"({p},{q}) dim {sdim}: dense λ = {lambda.Real:F6}{lambda.Imaginary:+0.000000;-0.000000}i");
                }
                bestOffset = Math.Min(bestOffset, (lambda - fold).Magnitude);

                var img = new Complex[tdim];
                ApplyLadder(n, p, q, phi: false, v, img);
                double nrm = Math.Sqrt(img.Sum(Abs2));
                if (Interior(n, p + 1, q - 1))
                {
                    var lt = BuildCsr(n, p + 1, q - 1, j);
                    var ltImg = new Complex[tdim];
                    CsrOps.Multiply(lt, img, ltImg);
                    double res = 0.0;
                    for (int i = 0; i < tdim; i++) res += Abs2(ltImg[i] - lambda * img[i]);
                    res = Math.Sqrt(res) / nrm;
                    rungs.Add(new SparseRung(p, q, lambda, nrm, res, sdim, tdim, sparse));
                }
                else
                {
                    terminal = nrm;
                    termSrc = sdim; termTgt = tdim;
                    // Survivors negative control: the terminal source is small at every walked N
                    // (324 at N=9), so the full dense eigenbasis stays affordable.
                    var a = new Complex[(long)sdim * sdim];
                    for (int r = 0; r < sdim; r++)
                        for (int k2 = ls.RowPtr[r]; k2 < ls.RowPtr[r + 1]; k2++)
                            a[r + (long)ls.ColIdx[k2] * sdim] = ls.Values[k2];
                    var (_, vrAll) = MklDirect.EigenvaluesAndVectorsDirectRaw(a, sdim);
                    var colv = new Complex[sdim];
                    var outv = new Complex[tdim];
                    for (int c = 0; c < sdim; c++)
                    {
                        double cn = 0.0;
                        for (int i = 0; i < sdim; i++) { colv[i] = vrAll[c * sdim + i]; cn += Abs2(colv[i]); }
                        cn = Math.Sqrt(cn);
                        ApplyLadder(n, p, q, phi: false, colv, outv);
                        if (Math.Sqrt(outv.Sum(Abs2)) / cn > 0.1) survivors++;
                    }
                }
            }
            pending.Add((total, rungs, terminal, survivors, termSrc, termTgt));
        }

        // FoldTolerance is INFORMATIONAL here, not a derived acceptance bound: unlike the dense walk,
        // which derives its tolerance from whole-block spectra, this walk only sees the selected fold
        // eigenvalues, so 10·min over them is circular as a gate (it always admits at least one rung).
        // The gate's actual bound lives in the test: 1e-2, backed by the measured distance of the
        // nearest UNRELATED eigenvalue (4.6e-2 at the 3024-dim N=9 block) against the pair members at
        // ≤ ~1.5e-3 from the recorded fold.
        double tol = Math.Max(bestOffset * 10, 1e-9);
        foreach (var (total, rungs, terminal, survivors, termSrc, termTgt) in pending)
        {
            double ratio = rungs.Count > 0 && !double.IsNaN(terminal)
                ? terminal / rungs.Max(r => r.Norm)
                : double.NaN;
            chainResults.Add(new SparseChain(total, rungs, terminal, ratio, survivors, termSrc, termTgt));
        }
        return new SparseWalkResult(n, qStar, lamA, fold, worstControl, tol, cg, chainResults);
    }
}
