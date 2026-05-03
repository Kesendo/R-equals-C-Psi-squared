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

        // 3c−2 modes is small at c=2..4 (4..10); serial Q loop is fine.
        for (int iQ = 0; iQ < qArr.Length; iQ++)
            ScanAtQ(iQ, qArr[iQ], tArr, kByBond, tAtPeak);

        return new KCurve(Effective.Block, qArr, tArr, kByBond, tAtPeak);
    }

    private double[] DefaultTGrid(int points = 21) =>
        ResonanceScan.DefaultTGrid(Effective.Block.GammaZero, points);

    private void ScanAtQ(int iQ, double q, double[] tGrid, double[,] kByBond, double[,] tAtPeak)
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

        for (int iT = 0; iT < tGrid.Length; iT++)
        {
            double t = tGrid[iT];
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
            var rhoT = R * weighted;

            for (int b = 0; b < numBonds; b++)
            {
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
                double kAbs = Math.Abs(2.0 * inner.Real);
                if (kAbs > kByBond[b, iQ])
                {
                    kByBond[b, iQ] = kAbs;
                    tAtPeak[b, iQ] = t;
                }
            }
        }
    }
}
