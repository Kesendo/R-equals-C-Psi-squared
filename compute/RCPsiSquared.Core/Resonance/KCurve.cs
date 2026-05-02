using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Resonance;

/// <summary>Per-bond K(Q) = max_t |∂S/∂J_b| curve over a Q grid. Carries the full per-bond
/// data and exposes bond-class averaging and peak/HWHM extraction.</summary>
public sealed record KCurve(
    CoherenceBlock Block,
    IReadOnlyList<double> QGrid,
    IReadOnlyList<double> TGrid,
    double[,] KByBond,
    double[,] TAtPeakByBond)
{
    public int NumBonds => KByBond.GetLength(0);

    public double[] BondClassAverage(BondClass cls)
    {
        var bonds = cls.BondsOf(NumBonds);
        int nQ = QGrid.Count;
        var avg = new double[nQ];
        if (bonds.Count == 0) return avg;
        for (int q = 0; q < nQ; q++)
        {
            double s = 0;
            foreach (int b in bonds) s += KByBond[b, q];
            avg[q] = s / bonds.Count;
        }
        return avg;
    }

    public double[] BondCurve(int bond)
    {
        int nQ = QGrid.Count;
        var c = new double[nQ];
        for (int q = 0; q < nQ; q++) c[q] = KByBond[bond, q];
        return c;
    }

    public PeakResult Peak(BondClass cls) => Peak(BondClassAverage(cls));
    public PeakResult PeakAtBond(int bond) => Peak(BondCurve(bond));

    public PeakResult Peak(double[] curve)
    {
        var info = ParabolicPeakFinder.Find(QGrid, curve);
        return new PeakResult(info.QPeak, info.KMax, info.HwhmLeft, info.HwhmRight);
    }
}
