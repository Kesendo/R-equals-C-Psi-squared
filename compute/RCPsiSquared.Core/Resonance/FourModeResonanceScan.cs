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

        // 4-mode is so small that parallelism per Q is overkill; the inner loop dominates anyway.
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

        const int dim = 4;
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
