using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Numerics;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Resonance;

/// <summary>F86 K_CC_pr Q-scan computed in the 4-mode effective basis (<see cref="FourModeEffective"/>).
///
/// <para>Same Duhamel + spectral-contour formulation as <see cref="ResonanceScan"/>, but on
/// 4×4 projected matrices instead of the full block-L. Output is a <see cref="KCurve"/>
/// compatible with the full scan, so direct shape comparison is straightforward.</para>
///
/// <para>The comparison answers the F86 open question: does the 4-mode minimal effective
/// model reproduce the universal shape K_class(x) = f_class(Q/Q_EP)? If yes, the closed
/// form for HWHM_left/Q_peak follows from the 4×4 eigenstructure analytically. If no, more
/// modes are needed — see PROOF_F86_QPEAK.md "What's missing for full Tier 1" item 1.</para>
/// </summary>
public sealed class FourModeResonanceScan
{
    public FourModeEffective Effective { get; }

    public FourModeResonanceScan(CoherenceBlock block) : this(FourModeEffective.Build(block)) { }

    public FourModeResonanceScan(FourModeEffective effective)
    {
        MathNetSetup.EnsureInitialized();
        Effective = effective;
    }

    public KCurve ComputeKCurve(IReadOnlyList<double>? qGrid = null, IReadOnlyList<double>? tGrid = null)
    {
        var qArr = (qGrid ?? ResonanceScan.DefaultQGrid()).ToArray();
        var tArr = (tGrid ?? DefaultTGrid()).ToArray();

        int numBonds = Effective.MhPerBondEff.Count;
        var kByBond = new double[numBonds, qArr.Length];
        var tAtPeak = new double[numBonds, qArr.Length];
        var atWindowEdge = new bool[numBonds, qArr.Length];

        // 4-mode is so small that parallelism per Q is overkill; the inner loop dominates anyway.
        for (int iQ = 0; iQ < qArr.Length; iQ++)
            ScanAtQ(iQ, qArr[iQ], tArr, kByBond, tAtPeak, atWindowEdge);

        return new KCurve(Effective.Block, qArr, tArr, kByBond, tAtPeak, atWindowEdge);
    }

    /// <summary>Points of this scan's default seed t grid: [0.2, 6]·t_peak at step 0.05·t_peak.</summary>
    public const int DefaultTPoints = 117;

    /// <summary>The four-mode model keeps the slow modes and drops the faster ones the full block carries,
    /// so its maximiser t*(Q) runs later than the full block's: at c = 2, γ₀ = 0.05 (t_peak = 5) the
    /// endpoint bond reaches t* = 16.4 (N = 5), 17.5 (N = 6), 23.8 (N = 7) and 18.5 (N = 8) at Q = 0.2,
    /// past the full block's window end 2.4·t_peak. This scan therefore seeds on [0.2, 6]·t_peak.</summary>
    public double[] DefaultTGrid(int points = DefaultTPoints)
    {
        double tPeak = EpAlgebra.TPeak(Effective.Block.GammaZero);
        return ResonanceScan.LinearQGrid(0.2 * tPeak, 6.0 * tPeak, points);
    }

    private void ScanAtQ(int iQ, double q, double[] tGrid, double[,] kByBond, double[,] tAtPeak, bool[,] atWindowEdge)
    {
        var L = Effective.LEffAtQ(q);
        var evd = L.Evd();
        var R = evd.EigenVectors;
        var Rinv = R.Inverse();
        var evals = evd.EigenValues.ToArray();

        var c0 = Rinv * Effective.ProbeEff;
        int numBonds = Effective.MhPerBondEff.Count;
        var xB = new ComplexMatrix[numBonds];
        for (int b = 0; b < numBonds; b++)
            xB[b] = Rinv * Effective.MhPerBondEff[b] * R;

        const int dim = 4;
        var expLam = new Complex[dim];
        var iMat = new Complex[dim, dim];
        var fbC0 = new Complex[dim];

        MathNet.Numerics.LinearAlgebra.Vector<Complex> rhoT = null!;

        // Fills expLam, iMat and rhoT for one t; KAbsAtBond then reads |K_b(Q, t)| off them.
        void Prepare(double t)
        {
            for (int i = 0; i < dim; i++) expLam[i] = Complex.Exp(evals[i] * t);

            for (int r = 0; r < dim; r++)
                for (int c = 0; c < dim; c++)
                {
                    Complex diff = evals[c] - evals[r];
                    iMat[r, c] = diff.Magnitude > 1e-10
                        ? (expLam[c] - expLam[r]) / diff
                        : t * expLam[r];
                }

            // ρ(t) = R · (expLam ⊙ c0)
            var weighted = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(dim);
            for (int i = 0; i < dim; i++) weighted[i] = expLam[i] * c0[i];
            rhoT = R * weighted;
        }

        double KAbsAtBond(int b)
        {
            // fbC0[r] = sum_c X_b[r, c] · iMat[r, c] · c0[c]
            for (int r = 0; r < dim; r++)
            {
                Complex s = Complex.Zero;
                for (int c = 0; c < dim; c++) s += xB[b][r, c] * iMat[r, c] * c0[c];
                fbC0[r] = s;
            }
            var fbVec = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(fbC0);
            var drho = R * fbVec;
            var sDrho = Effective.SKernelEff * drho;
            Complex inner = Complex.Zero;
            for (int i = 0; i < dim; i++) inner += Complex.Conjugate(rhoT[i]) * sDrho[i];
            return Math.Abs(2.0 * inner.Real);
        }

        // Seed: the grid maximum per bond; then the envelope, refined in continuous t.
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
}
