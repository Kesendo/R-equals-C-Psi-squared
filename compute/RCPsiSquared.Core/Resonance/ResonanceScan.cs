using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Probes;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Resonance;

/// <summary>Per-bond K_CC_pr = max_t |∂S/∂J_b|(Q) scanner.
///
/// For uniform J at each Q: build L(Q) = D + J·M_H_total, diagonalize, evaluate ρ(t)
/// and ∂ρ/∂J_b via Duhamel, observable K_b(Q, t) = 2·Re(⟨ρ(t)| S_kernel | ∂ρ/∂J_b⟩).
/// Peak over t per (bond, Q): the t-grid seeds the maximum and <see cref="EnvelopeRefinement"/>
/// refines it in continuous t, so the curve is the envelope max_t, not a grid reading.
///
/// Performance: Q-grid parallelised via Parallel.For; the heavy O(dim³) operations
/// (L assembly, eigendecomposition, inverse, X_b = R⁻¹ · M_H_b · R) go through MathNet
/// matrix ops which dispatch to native MKL once <see cref="MathNetSetup"/> is initialised.
/// Per-thread scratch is reused via <see cref="ThreadLocal{T}"/>.
/// </summary>
public sealed class ResonanceScan
{
    public CoherenceBlock Block { get; }
    public ComplexVector Probe { get; }
    public ComplexMatrix SKernel { get; }

    public ResonanceScan(CoherenceBlock block)
    {
        MathNetSetup.EnsureInitialized();
        Block = block;
        Probe = DickeBlockProbe.Build(block);
        SKernel = SpatialSumKernel.Build(block);
    }

    /// <summary>Default Q grid: dQ = 0.025, range [0.20, 4.00] (153 points).</summary>
    public static double[] DefaultQGrid() => LinearQGrid(0.20, 4.00, 153);

    /// <summary>Reduced-resolution Q grid for fast inspection: <paramref name="points"/>
    /// linearly spaced between <paramref name="lo"/> and <paramref name="hi"/>.</summary>
    public static double[] LinearQGrid(double lo, double hi, int points)
    {
        if (points < 2) throw new ArgumentOutOfRangeException(nameof(points), $"need ≥2 points; got {points}");
        var g = new double[points];
        for (int i = 0; i < points; i++) g[i] = lo + (hi - lo) * i / (points - 1);
        return g;
    }

    /// <summary>Number of points of the default seed t grid (step 0.05·t_peak).</summary>
    public const int DefaultTPoints = 45;

    /// <summary>Default seed t grid: <see cref="DefaultTPoints"/> points spanning [0.2, 2.4]·t_peak
    /// with t_peak = 1/(4γ₀). It seeds the maximum and bounds the window; the value is refined in
    /// continuous t (<see cref="EnvelopeRefinement"/>).</summary>
    public double[] DefaultTGrid(int points = DefaultTPoints) => DefaultTGrid(Block.GammaZero, points);

    /// <summary>Default seed t grid for any γ₀, shared between full-block and effective-mode scanners.
    ///
    /// <para>The window has to contain the envelope's maximiser t*(Q) wherever a curve is read,
    /// the half-maximum flanks included, or the read is a window edge rather than the envelope.
    /// Measured at c = 2, N = 6, γ₀ = 0.05 (t_peak = 5): t*(Q) approaches 2·t_peak = 1/(2γ₀) as Q → 0
    /// (9.88 at Q = 0.2, 9.99 at Q = 0.05) and falls to about 0.36·t_peak at large Q (1.81 at Q = 5), with branch
    /// jumps where K(t) is multimodal. [0.2, 2.4]·t_peak holds all of these at N = 6. It does not hold
    /// every branch at larger N: the third bond from each end carries a late branch on the far right
    /// flank (bond 2): t* = 11.94 at Q = 3.85 for N = 10, just inside the window's end at 12, and
    /// 12.68 at Q = 3.65 for N = 11 and 12.93 for N = 12, past it. Those cells lie far above that
    /// bond's Q_peak (1.53) and do not reach its
    /// Q_peak or its left half-maximum; a cell whose grid winner sits on a window end is flagged in
    /// <see cref="KCurve.TStarAtWindowEdge"/>.</para></summary>
    public static double[] DefaultTGrid(double gammaZero, int points = DefaultTPoints)
    {
        double tPeak = EpAlgebra.TPeak(gammaZero);
        return LinearQGrid(0.2 * tPeak, 2.4 * tPeak, points);
    }

    /// <summary>The K-curve over the Q grid, one Q point per parallel work item.
    /// <paramref name="maxDegreeOfParallelism"/> caps the workers (−1, the default, lets the
    /// runtime choose); each worker holds its own dim × dim workspace, so at large blocks
    /// (c = 4, N = 8, dim 3920) the cap is what keeps the scan inside memory.</summary>
    public KCurve ComputeKCurve(IReadOnlyList<double>? qGrid = null, IReadOnlyList<double>? tGrid = null,
        int maxDegreeOfParallelism = -1)
    {
        var qArr = (qGrid ?? DefaultQGrid()).ToArray();
        var tArr = (tGrid ?? DefaultTGrid()).ToArray();

        var ctx = new ScanContext(Block, Probe, SKernel, tArr);
        int numBonds = ctx.Decomp.NumBonds;
        var kByBond = new double[numBonds, qArr.Length];
        var tAtPeak = new double[numBonds, qArr.Length];
        var atWindowEdge = new bool[numBonds, qArr.Length];

        using var workPool = new ThreadLocal<ScanWorkspace>(() => new ScanWorkspace(ctx.Dim, numBonds), trackAllValues: false);
        var options = new ParallelOptions { MaxDegreeOfParallelism = maxDegreeOfParallelism };
        Parallel.For(0, qArr.Length, options, iQ =>
        {
            ScanAtQ(in ctx, qArr[iQ], iQ, workPool.Value!, kByBond, tAtPeak, atWindowEdge);
        });

        return new KCurve(Block, qArr, tArr, kByBond, tAtPeak, atWindowEdge);
    }

    private static void ScanAtQ(in ScanContext ctx, double q, int iQ, ScanWorkspace work,
        double[,] kByBond, double[,] tAtPeak, bool[,] atWindowEdge)
    {
        int dim = ctx.Dim;
        int numBonds = ctx.Decomp.NumBonds;
        double j = q * ctx.GammaZero;

        // L(Q) = D + j · MhTotal (MKL via MathNet matrix ops)
        var L = ctx.Decomp.D + (Complex)j * ctx.Decomp.MhTotal;
        var evd = L.Evd();
        var R = evd.EigenVectors;
        var Rinv = R.Inverse();

        // c0 = R⁻¹ · probe (MKL matvec)
        var c0 = Rinv * ctx.Probe;

        // Per-bond X_b = R⁻¹ · M_H_per_bond[b] · R (MKL gemm) — the O(dim³) work that previously
        // bypassed MKL. Extract to Complex[,] for the element-wise inner product below.
        var xBArr = work.XbArrays;
        for (int b = 0; b < numBonds; b++)
        {
            var xB = Rinv * ctx.Decomp.MhPerBond[b] * R;
            xBArr[b] = xB.ToArray();
        }

        var rArr = R.ToArray();
        var c0Arr = c0.ToArray();
        var evals = evd.EigenValues.ToArray();
        var iMat = work.IMat;
        var expLam = work.ExpLam;
        var rhoT = work.RhoT;
        var fbC0 = work.FbC0;
        var sKernel = ctx.SKernel;   // an `in` parameter cannot be captured by the local functions
        var tGrid = ctx.TGrid;

        // Fills expLam, iMat and rhoT for one t; KAbsAtBond then reads |K_b(Q, t)| off them.
        void Prepare(double t)
        {
            for (int i = 0; i < dim; i++) expLam[i] = Complex.Exp(evals[i] * t);

            // I_mat[r, c] = (e_c - e_r) / (λ_c - λ_r) or t·e_r at degeneracy
            for (int r = 0; r < dim; r++)
            {
                for (int c = 0; c < dim; c++)
                {
                    Complex diff = evals[c] - evals[r];
                    iMat[r, c] = diff.Magnitude > 1e-10
                        ? (expLam[c] - expLam[r]) / diff
                        : t * expLam[r];
                }
            }

            // ρ(t) = R · (expLam ⊙ c0)
            for (int i = 0; i < dim; i++)
            {
                Complex s = Complex.Zero;
                for (int kk = 0; kk < dim; kk++) s += rArr[i, kk] * (expLam[kk] * c0Arr[kk]);
                rhoT[i] = s;
            }
        }

        double KAbsAtBond(int b)
        {
            // fbC0[r] = sum_c X_b[r, c] · iMat[r, c] · c0[c]   (no `weighted` matrix; fold c0 into the inner sum)
            var xB = xBArr[b];
            for (int r = 0; r < dim; r++)
            {
                Complex s = Complex.Zero;
                for (int c = 0; c < dim; c++) s += xB[r, c] * iMat[r, c] * c0Arr[c];
                fbC0[r] = s;
            }
            // drho = R · fbC0  (MKL matvec via Vector<Complex>)
            var fbC0Vec = ComplexVector.Build.Dense(fbC0);
            var drho = R * fbC0Vec;
            // sDrho = SKernel · drho (MKL matvec)
            var sDrho = sKernel * drho;
            // K = 2 · Re(⟨rhoT | sDrho⟩)
            Complex inner = Complex.Zero;
            for (int i = 0; i < dim; i++) inner += Complex.Conjugate(rhoT[i]) * sDrho[i];
            return Math.Abs(2.0 * inner.Real);
        }

        // Seed: the grid maximum per bond.
        var winIdx = new int[numBonds];
        var winK = new double[numBonds];
        for (int iT = 0; iT < tGrid.Length; iT++)
        {
            Prepare(tGrid[iT]);
            for (int b = 0; b < numBonds; b++)
            {
                double kAbs = KAbsAtBond(b);
                if (kAbs > winK[b]) { winK[b] = kAbs; winIdx[b] = iT; }
            }
        }

        // Envelope: refine each bond's maximum in continuous t around its seed.
        for (int b = 0; b < numBonds; b++)
        {
            int bond = b;
            var (kMax, tStar, atEdge) = EnvelopeRefinement.Refine(tGrid, winIdx[b], winK[b],
                t => { Prepare(t); return KAbsAtBond(bond); });
            kByBond[b, iQ] = kMax;
            tAtPeak[b, iQ] = tStar;
            atWindowEdge[b, iQ] = atEdge;
        }
    }

    private readonly struct ScanContext
    {
        public readonly CoherenceBlock Block;
        public readonly BlockLDecomposition Decomp;
        public readonly ComplexVector Probe;
        public readonly ComplexMatrix SKernel;
        public readonly double[] TGrid;
        public readonly double GammaZero;
        public readonly int Dim;

        public ScanContext(CoherenceBlock block, ComplexVector probe, ComplexMatrix sKernel, double[] tGrid)
        {
            Block = block;
            Decomp = block.Decomposition;
            Probe = probe;
            SKernel = sKernel;
            TGrid = tGrid;
            GammaZero = block.GammaZero;
            Dim = block.Basis.MTotal;
        }
    }

    private sealed class ScanWorkspace
    {
        public readonly Complex[,] IMat;
        public readonly Complex[] ExpLam;
        public readonly Complex[] RhoT;
        public readonly Complex[] FbC0;
        public readonly Complex[][,] XbArrays;

        public ScanWorkspace(int dim, int numBonds)
        {
            IMat = new Complex[dim, dim];
            ExpLam = new Complex[dim];
            RhoT = new Complex[dim];
            FbC0 = new Complex[dim];
            XbArrays = new Complex[numBonds][,];
        }
    }
}
