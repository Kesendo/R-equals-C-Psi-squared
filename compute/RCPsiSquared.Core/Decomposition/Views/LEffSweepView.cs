using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Core.Decomposition.Views;

/// <summary>F86 4-mode L_eff(Q) sweep over a Q grid: the 3D / 4D structure that the Object
/// Manager exports.
///
/// <para>Dimensions: <c>Q × {L_eff matrix 4×4} × {eigenvalues 4-vec, eigenvectors 4×4}</c>.
/// Every Q-slice is a <see cref="LEffSpectralAtQ"/>, itself a computational unit (lazy EVD,
/// max-imaginary-part as its scalar). The composition with
/// <c>FourModeEffective.MhPerBondViews</c> at the parent level pushes this into 4D
/// (bond × Q × 4 × 4).</para>
///
/// <para>Also computes a Q-vs-max|Im(λ)| curve as a curve payload. The curve is non-zero from the
/// smallest Q (the Hamiltonian part oscillates on its own) and grows with Q; its slope can change where
/// a pair leaves the real axis. It does not peak at an exceptional point, so its maximum over a grid
/// sits at the grid's upper end (Q = 4.0 on the default grid at N = 5, 6, 7) and locates no EP.</para>
/// </summary>
public sealed class LEffSweepView : IInspectable
{
    public FourModeEffective Effective { get; }
    public IReadOnlyList<double> QGrid { get; }

    private readonly Lazy<IReadOnlyList<LEffSpectralAtQ>> _snapshots;
    private readonly Lazy<double[]> _maxImCurve;
    private readonly Lazy<double> _qAtPeakIm;

    public IReadOnlyList<LEffSpectralAtQ> Snapshots => _snapshots.Value;

    /// <summary>max |Im(λ)| as a function of Q over the grid.</summary>
    public IReadOnlyList<double> MaxImaginaryCurve => _maxImCurve.Value;

    /// <summary>The Q value where max |Im(λ)| is largest in the grid; the curve grows with Q, so this is the
    /// grid's upper end in every sweep read so far, not an exceptional point.</summary>
    public double QAtMaxImaginary => _qAtPeakIm.Value;

    public LEffSweepView(FourModeEffective effective, IReadOnlyList<double> qGrid)
    {
        Effective = effective;
        QGrid = qGrid;

        _snapshots = new Lazy<IReadOnlyList<LEffSpectralAtQ>>(() =>
        {
            var arr = new LEffSpectralAtQ[qGrid.Count];
            for (int i = 0; i < qGrid.Count; i++)
                arr[i] = new LEffSpectralAtQ(qGrid[i], effective.Block.GammaZero, effective.LEffAtQ(qGrid[i]));
            return arr;
        });
        _maxImCurve = new Lazy<double[]>(() =>
        {
            var c = new double[qGrid.Count];
            for (int i = 0; i < qGrid.Count; i++) c[i] = Snapshots[i].MaxImaginaryPart;
            return c;
        });
        _qAtPeakIm = new Lazy<double>(() =>
        {
            int iMax = 0;
            for (int i = 1; i < MaxImaginaryCurve.Count; i++)
                if (MaxImaginaryCurve[i] > MaxImaginaryCurve[iMax]) iMax = i;
            return qGrid[iMax];
        });
    }

    public string DisplayName => $"L_eff sweep over Q ({QGrid.Count} points)";
    public string Summary =>
        $"Q ∈ [{QGrid.Min():F3}, {QGrid.Max():F3}], max |Im(λ)| largest at Q = {QAtMaxImaginary:F4}{(QAtMaxImaginary == QGrid[^1] ? " (the grid's end, no EP)" : " (inside the grid)")}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("Q at peak |Im(λ)|", QAtMaxImaginary, "F4");
            yield return new InspectableNode("max |Im(λ)| curve",
                summary: $"max over Q = {MaxImaginaryCurve.Max():E3}",
                payload: new InspectablePayload.Curve("max|Im(λ)|(Q)", QGrid, MaxImaginaryCurve.ToArray(),
                    XLabel: "Q", YLabel: "max |Im(λ_eff)|"));
            yield return InspectableNode.Group("snapshots", Snapshots, Snapshots.Count);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
