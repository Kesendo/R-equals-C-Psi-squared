using System.Numerics;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The sideways-ladder chain walk of the <c>sideways_spin_ladder</c> open arc: the ONE
/// construction behind both the C# gate (<c>SidewaysSpinLadderGateTests</c>, Categories SIDEWAYS +
/// SLOW_SIDEWAYS) and the live witness (<see cref="SidewaysSpinLadderWitness"/>), so the gate's numbers
/// and the inspect-time numbers cannot drift apart. Port of <c>simulations/eta_ladder_chain.py</c>.
///
/// <para>What one <see cref="Run"/> computes, in the order that makes the later parts readable:
/// CONTROL — Φ (<see cref="SpectatorIntertwiner.BuildW"/>) and S⁺
/// (<see cref="SpectatorIntertwiner.BuildSPlus"/>) both intertwine L on every rung of the two chains
/// p+q̃ = N∓1, worst residual compared to 0.0 exactly (the Python gate measures exactly 0.0 here, and so
/// does this one; see <see cref="ExactResidual"/> for why the zero is order-independent). SHAPE — the interior blocks whose spectrum carries
/// the cross-fold partner −λ_A − 2N, against the two chain interiors. COUNT — fold + band sectors vs the
/// F125 orbit size 4N−8. LADDER — per chain, the transport norms ‖S⁺v‖ of the fold eigenvector against
/// the Clebsch-Gordan coefficients √(ℓ(ℓ+1) − m(m+1)), ℓ = (N−3)/2, plus the terminal step into the
/// boundary block as a RATIO to the interior norms (the vector comes from an eigensolver, so the exact
/// zero of the highest-weight death S⁺|ℓ,ℓ⟩ = 0 is not measurable, its floor ratio is) and the survivors
/// negative control (the terminal map is not the zero map).</para>
///
/// <para><b>Conventions</b> (identical to <c>SpectatorIntertwinerGateTests</c>): H = XYChain(N, 2q*) so
/// the hop amplitude is J = 2q*, γ = 1 per site, basis = <see cref="BlockBasis.PopcountStates"/>
/// ascending, flat = pIdx·Mq + qIdx, blocks via <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/>,
/// eigensolver = LAPACK zgeev through <see cref="MklDirect"/>. The recorded λ_A is a 4-decimal value, so
/// the acceptance tolerance is DERIVED from the run (ten times the smallest offset actually found)
/// rather than set.</para></summary>
public static class SidewaysSpinLadderChain
{
    /// <summary>One interior rung of a chain: the fold eigenvector of (P,Q) transported to
    /// (P+1,Q−1) at norm <paramref name="Norm"/>, staying an eigenvector at the same λ up to
    /// <paramref name="EigResidual"/> (relative). <paramref name="TargetDim"/> is the target block's
    /// dimension, the scale of the eigensolver error model the gate normalizes the residual by.</summary>
    public readonly record struct Rung(int P, int Q, double Norm, double EigResidual, int TargetDim);

    /// <summary>One chain p+q̃ = Total: the interior rungs, the terminal norm into the boundary
    /// block, its ratio to the largest interior norm, and the survivors negative control
    /// (how many of the terminal source block's eigenvectors transport at norm &gt; 0.1;
    /// must be ≥ TerminalTargetDim for the death of ONE vector to mean anything).</summary>
    public readonly record struct Chain(
        int Total, IReadOnlyList<Rung> Rungs, double Terminal, double TerminalRatio,
        int Survivors, int TerminalSourceDim, int TerminalTargetDim);

    /// <summary>Everything one run at (N, q*, λ_A) produces; the gate asserts on these fields, the
    /// witness renders them.</summary>
    public sealed record LadderRun(
        int N, double QStar, double LambdaA, double Fold,
        double WorstControlResidual, double Tolerance,
        IReadOnlyList<(int P, int Q)> FoldSet, IReadOnlyList<(int P, int Q)> BandSet,
        IReadOnlyList<(int P, int Q)> PredictedSet,
        IReadOnlyList<double> CgPredicted, IReadOnlyList<Chain> Chains);

    /// <summary>The recorded seed used per N: the same census entries the Python gate hard-codes
    /// (N=5 locus 1 and the N=7 q*≈1.515 locus, both R-parity +1). Verified against
    /// <see cref="RealDefectiveSeeds.All"/> by the gate's seed test.</summary>
    public static (double QStar, double LambdaA) Seed(int n) => n switch
    {
        5 => (0.620878, -4.6189),
        7 => (1.514833, -4.8846),
        _ => throw new ArgumentOutOfRangeException(nameof(n),
            $"no seed wired for N={n} (RealDefectiveSeeds lists odd N; the dense chain walk is seeded " +
            "at N=5 and N=7, the N=9 seed lives in SidewaysSpinLadderSparse.Seed)"),
    };

    /// <summary>The Clebsch-Gordan transport norms √(ℓ(ℓ+1) − m(m+1)) for ℓ = (N−3)/2,
    /// m = −ℓ … ℓ−1 dropped one short (the interior rungs; length N−3).</summary>
    public static IReadOnlyList<double> CgCoefficients(int n)
    {
        double spin = (n - 3) / 2.0;
        var cg = new List<double>();
        for (double m = -spin; m < spin - 0.5; m += 1.0)
            cg.Add(Math.Sqrt(spin * (spin + 1) - m * (m + 1)));
        return cg;
    }

    /// <summary>The (p,q̃) joint-popcount block of L at real q: H = XYChain(n, 2q), γ = 1 per site,
    /// flat = pIdx·Mq + qIdx, exactly as <c>SpectatorIntertwinerGateTests.BuildBlockN</c>.</summary>
    public static ComplexMatrix BuildBlock(int n, int p, int qt, double q)
        => BuildBlockWithH(PauliHamiltonian.XYChain(n, 2.0 * q).ToMatrix(), n, p, qt);

    /// <summary>Same block, from a prebuilt H (one Run touches ~200 blocks; H is built once).</summary>
    private static ComplexMatrix BuildBlockWithH(ComplexMatrix h, int n, int p, int qt)
    {
        int d = 1 << n;
        var statesP = BlockBasis.PopcountStates(n, p);
        var statesQ = BlockBasis.PopcountStates(n, qt);
        var gamma = Enumerable.Repeat(1.0, n).ToList();
        var flat = new int[statesP.Count * statesQ.Count];
        for (int i = 0; i < statesP.Count; i++)
            for (int j = 0; j < statesQ.Count; j++)
                flat[i * statesQ.Count + j] = (int)(statesP[i] * d + statesQ[j]);
        return PerBlockLiouvillianBuilder.BuildBlockZ(h, gamma, flat);
    }

    /// <summary>Max-entry magnitude of L_target·M − M·L_source (max-entry rather than Python's Frobenius:
    /// identical as an exact-zero test, different as a diagnostic magnitude if it ever went nonzero). The
    /// comparison to 0.0 is the house case-1 rule: the operator identity cancels pairwise per entry, both
    /// sides of each entry draw the same multiset of addends (products of L entries with the ±1 ladder
    /// entries; γ = 1 real parts are integer multiples of 2, imaginary parts integer multiples of the one
    /// hop amplitude J), and the measured residual IS exactly 0.0, in numpy's zgemm and in this ordered
    /// walk alike, on every rung at N = 5 and 7. A nonzero here is a finding about the construction, not a
    /// tolerance case. Sparse column walk: every column of M has ≤ n nonzeros and every column of L has
    /// ≤ 2(n−1)+1, so the residual costs O(dim²·n), not O(dim³).</summary>
    public static double ExactResidual(Complex[,] lt, Complex[,] m, Complex[,] ls)
    {
        int tdim = lt.GetLength(0), sdim = ls.GetLength(0);
        double worst = 0.0;
        var colA = new Complex[tdim];
        var colB = new Complex[tdim];
        for (int j = 0; j < sdim; j++)
        {
            Array.Clear(colA);
            Array.Clear(colB);
            for (int k = 0; k < tdim; k++)          // (Lt·M)[:,j] = Σ_k M[k,j]·Lt[:,k], ascending k
            {
                var mkj = m[k, j];
                if (mkj == Complex.Zero) continue;
                for (int i = 0; i < tdim; i++)
                    if (lt[i, k] != Complex.Zero)
                        colA[i] += lt[i, k] * mkj;
            }
            for (int k = 0; k < sdim; k++)          // (M·Ls)[:,j] = Σ_k Ls[k,j]·M[:,k], ascending k
            {
                var lkj = ls[k, j];
                if (lkj == Complex.Zero) continue;
                for (int i = 0; i < tdim; i++)
                    if (m[i, k] != Complex.Zero)
                        colB[i] += m[i, k] * lkj;
            }
            for (int i = 0; i < tdim; i++)
                worst = Math.Max(worst, (colA[i] - colB[i]).Magnitude);
        }
        return worst;
    }

    private static int Binom(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        long r = 1;
        for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
        return (int)r;
    }

    private static bool Interior(int n, int p, int q) => p >= 1 && p <= n - 1 && q >= 1 && q <= n - 1;

    private static double Abs2(Complex x) => x.Real * x.Real + x.Imaginary * x.Imaginary;

    private static Complex[] ToColMajor(Complex[,] m)
    {
        int rows = m.GetLength(0), cols = m.GetLength(1);
        var a = new Complex[(long)rows * cols];
        for (int j = 0; j < cols; j++)
            for (int i = 0; i < rows; i++)
                a[i + (long)j * rows] = m[i, j];
        return a;
    }

    /// <summary>The full chain walk at the recorded seed for <paramref name="n"/> (see
    /// <see cref="Seed"/>), or at an explicit (q*, λ_A).</summary>
    public static LadderRun Run(int n)
    {
        var (qStar, lamA) = Seed(n);
        return Run(n, qStar, lamA);
    }

    /// <summary>The full chain walk at an explicit (q*, λ_A): CONTROL residuals, the per-block spectra
    /// with the derived tolerance, the fold/band sector sets, and both chains' transport norms with
    /// terminal ratio and survivors control. The gate asserts on the returned fields; the witness
    /// renders them.</summary>
    public static LadderRun Run(int n, double qStar, double lamA)
    {
        double fold = -lamA - 2.0 * n;
        int[] chains = { n - 1, n + 1 };
        var h = PauliHamiltonian.XYChain(n, 2.0 * qStar).ToMatrix();   // built once, ~200 blocks reuse it

        // CONTROL: both ladders intertwine exactly, on every rung of both chains.
        double worstControl = 0.0;
        foreach (bool phi in new[] { true, false })
        {
            int dq = phi ? 1 : -1;
            foreach (int total in chains)
                for (int p = 0; p <= n; p++)
                {
                    int q = total - p;
                    if (q < 0 || q > n || p + 1 > n || q + dq < 0 || q + dq > n) continue;
                    int sdim = Binom(n, p) * Binom(n, q);
                    int tdim = Binom(n, p + 1) * Binom(n, q + dq);
                    if (sdim == 0 || tdim == 0) continue;
                    var ls = BuildBlockWithH(h, n, p, q).ToArray();
                    var lt = BuildBlockWithH(h, n, p + 1, q + dq).ToArray();
                    var m = (phi ? SpectatorIntertwiner.BuildW(n, p, q)
                                 : SpectatorIntertwiner.BuildSPlus(n, p, q)).ToArray();
                    worstControl = Math.Max(worstControl, ExactResidual(lt, m, ls));
                }
        }

        // Spectra once per block, and a tolerance derived from the run (λ_A is 4-decimal).
        var spec = new Dictionary<(int P, int Q), Complex[]>();
        double bestOffset = double.PositiveInfinity;
        for (int p = 0; p <= n; p++)
            for (int q = 0; q <= n; q++)
            {
                int dim = Binom(n, p) * Binom(n, q);
                if (dim == 0) continue;
                var eigs = MklDirect.EigenvaluesOnlyRaw(ToColMajor(BuildBlockWithH(h, n, p, q).ToArray()), dim);
                spec[(p, q)] = eigs;
                bestOffset = Math.Min(bestOffset, eigs.Min(e => (e - fold).Magnitude));
                bestOffset = Math.Min(bestOffset, eigs.Min(e => (e - lamA).Magnitude));
            }
        double tol = Math.Max(bestOffset * 10, 1e-9);

        var foldSet = spec.Where(kv => Interior(n, kv.Key.P, kv.Key.Q)
                                       && kv.Value.Min(e => (e - fold).Magnitude) < tol)
                          .Select(kv => kv.Key).OrderBy(k => k).ToList();
        var bandSet = spec.Where(kv => Interior(n, kv.Key.P, kv.Key.Q)
                                       && kv.Value.Min(e => (e - lamA).Magnitude) < tol)
                          .Select(kv => kv.Key).OrderBy(k => k).ToList();
        var predSet = chains.SelectMany(total => Enumerable.Range(0, n + 1)
                                .Where(p => total - p >= 0 && total - p <= n && Interior(n, p, total - p))
                                .Select(p => (P: p, Q: total - p)))
                            .OrderBy(k => k).ToList();

        // LADDER: transport norms per chain against the Clebsch-Gordan coefficients.
        var chainResults = new List<Chain>();
        foreach (int total in chains)
        {
            var walk = Enumerable.Range(0, n + 1).Where(p => total - p >= 0 && total - p <= n)
                                 .Select(p => (P: p, Q: total - p)).ToList();
            var rungs = new List<Rung>();
            double terminal = double.NaN;
            Complex[,]? termLs = null, termS = null;
            int termTgtDim = 0, termSrcDim = 0;
            foreach (var (p, q) in walk.Take(walk.Count - 1))
            {
                int sdim = Binom(n, p) * Binom(n, q);
                int tdim = Binom(n, p + 1) * Binom(n, q - 1);
                if (sdim == 0 || tdim == 0) continue;
                var lsArr = BuildBlockWithH(h, n, p, q).ToArray();
                var (w, vr) = MklDirect.EigenvaluesAndVectorsDirectRaw(ToColMajor(lsArr), sdim);
                int jBest = 0;
                for (int j = 1; j < sdim; j++)
                    if ((w[j] - fold).Magnitude < (w[jBest] - fold).Magnitude) jBest = j;
                if ((w[jBest] - fold).Magnitude > tol) continue;   // source does not carry the fold
                var v = new Complex[sdim];
                double vn = 0.0;
                for (int i = 0; i < sdim; i++) { v[i] = vr[jBest * sdim + i]; vn += Abs2(v[i]); }
                vn = Math.Sqrt(vn);
                for (int i = 0; i < sdim; i++) v[i] /= vn;
                var s = SpectatorIntertwiner.BuildSPlus(n, p, q).ToArray();
                var img = new Complex[tdim];
                for (int k = 0; k < sdim; k++)
                {
                    if (v[k] == Complex.Zero) continue;
                    for (int i = 0; i < tdim; i++)
                        if (s[i, k] != Complex.Zero)
                            img[i] += s[i, k] * v[k];
                }
                double nrm = Math.Sqrt(img.Sum(Abs2));
                if (Interior(n, p + 1, q - 1))
                {
                    var ltArr = BuildBlockWithH(h, n, p + 1, q - 1).ToArray();
                    double res = 0.0;
                    for (int i = 0; i < tdim; i++)
                    {
                        Complex acc = -w[jBest] * img[i];
                        for (int k = 0; k < tdim; k++) acc += ltArr[i, k] * img[k];
                        res += Abs2(acc);
                    }
                    res = Math.Sqrt(res) / nrm;
                    rungs.Add(new Rung(p, q, nrm, res, tdim));
                }
                else
                {
                    terminal = nrm;
                    termLs = lsArr; termS = s; termTgtDim = tdim; termSrcDim = sdim;
                }
            }

            // The terminal step is NOT an exact zero and cannot be one: the intertwining residual IS
            // exactly 0.0 (operator identity, cancels pairwise), but v comes from an eigensolver, so the
            // meaningful quantity is the RATIO to the interior steps. The terminal check is weak on its
            // own (most of the source block is in the kernel at the last rung); the survivors negative
            // control is what makes the death readable.
            int survivors = 0;
            if (termLs is not null && termS is not null)
            {
                var (_, vrAll) = MklDirect.EigenvaluesAndVectorsDirectRaw(ToColMajor(termLs), termSrcDim);
                for (int j = 0; j < termSrcDim; j++)
                {
                    double cn = 0.0;
                    for (int i = 0; i < termSrcDim; i++) cn += Abs2(vrAll[j * termSrcDim + i]);
                    cn = Math.Sqrt(cn);
                    double outn = 0.0;
                    for (int i = 0; i < termTgtDim; i++)
                    {
                        Complex acc = Complex.Zero;
                        for (int k = 0; k < termSrcDim; k++)
                            if (termS[i, k] != Complex.Zero)
                                acc += termS[i, k] * vrAll[j * termSrcDim + k];
                        outn += Abs2(acc);
                    }
                    if (Math.Sqrt(outn) / cn > 0.1) survivors++;
                }
            }
            double ratio = rungs.Count > 0 && !double.IsNaN(terminal)
                ? terminal / rungs.Max(r => r.Norm)
                : double.NaN;
            chainResults.Add(new Chain(total, rungs, terminal, ratio, survivors, termSrcDim, termTgtDim));
        }

        return new LadderRun(n, qStar, lamA, fold, worstControl, tol,
            foldSet, bandSet, predSet, CgCoefficients(n), chainResults);
    }
}
