using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Numerics;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Resonance;

/// <summary>F86 Item 2: K_CC_pr Q-scan in the multi-k effective basis (3c−2 modes).
/// Generalises <see cref="FourModeResonanceScan"/> to chromaticity c ≥ 2 — for c=2 these
/// give identical results.
///
/// <para>Comparison vs full block-L <see cref="ResonanceScan"/> answers F86 Open Item 2:
/// does the multi-k effective reproduce the K-curve at c≥3? If yes, the slowest-pair (k=1)
/// quartet structure extends naturally; if no, even more modes are needed.</para>
/// </summary>
public sealed class MultiKResonanceScan
{
    public MultiKEffective Effective { get; }

    public MultiKResonanceScan(CoherenceBlock block) : this(MultiKEffective.Build(block)) { }

    public MultiKResonanceScan(MultiKEffective effective)
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

        // 3c−2 modes is small at c=2..4 (4..10); serial Q loop is fine.
        for (int iQ = 0; iQ < qArr.Length; iQ++)
            ScanAtQ(iQ, qArr[iQ], tArr, kByBond, tAtPeak, atWindowEdge);

        return new KCurve(Effective.Block, qArr, tArr, kByBond, tAtPeak, atWindowEdge);
    }

    /// <summary>The seed t window of the effective-basis scans, [0.2, 6]·t_peak with
    /// <see cref="FourModeResonanceScan.DefaultTPoints"/> points: like the four-mode model the
    /// multi-k basis keeps only slow modes, so its t* runs past the full block's window end
    /// 2.4·t_peak at low Q (the four-mode values are in <see cref="FourModeResonanceScan"/>).</summary>
    private double[] DefaultTGrid(int points = FourModeResonanceScan.DefaultTPoints)
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

        int dim = Effective.Dim;
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
