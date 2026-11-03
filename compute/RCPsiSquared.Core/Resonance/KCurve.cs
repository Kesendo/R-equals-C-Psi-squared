using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Resonance;

/// <summary>Per-bond K(Q) = max_t |∂S/∂J_b| curve over a Q grid, the F86 envelope: <see cref="TGrid"/>
/// is the seed grid and window, <see cref="TAtPeakByBond"/> the refined maximiser t* (see
/// <see cref="EnvelopeRefinement"/>). <see cref="TStarAtWindowEdge"/> marks (bond, Q) cells whose
/// grid winner sat on a window end, where the maximiser may lie outside the window. Carries the
/// full per-bond data and exposes bond-class averaging and peak/HWHM extraction.</summary>
public sealed record KCurve(
    CoherenceBlock Block,
    IReadOnlyList<double> QGrid,
    IReadOnlyList<double> TGrid,
    double[,] KByBond,
    double[,] TAtPeakByBond,
    bool[,] TStarAtWindowEdge)
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

    /// <summary>The canonical peak of a bond class: the first local maximum from below
    /// (<see cref="ParabolicPeakFinder.FindFirstLocalMaximum"/>). <see cref="Peak(BondClass)"/>
    /// is the global maximum over this curve's Q window, which moves when the window reaches a
    /// bond that rises again at large Q.</summary>
    public PeakResult CanonicalPeak(BondClass cls) => CanonicalPeak(BondClassAverage(cls));
    public PeakResult CanonicalPeakAtBond(int bond) => CanonicalPeak(BondCurve(bond));

    public PeakResult CanonicalPeak(double[] curve)
    {
        var info = ParabolicPeakFinder.FindFirstLocalMaximum(QGrid, curve);
        return new PeakResult(info.QPeak, info.KMax, info.HwhmLeft, info.HwhmRight);
    }
}
