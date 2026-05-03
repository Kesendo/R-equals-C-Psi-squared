using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using RCPsiSquared.Core.Inspection;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Decomposition.Views;

/// <summary>F86 4-mode L_eff(Q) snapshot at one Q. Computes the matrix L_eff = D_eff +
/// (Q·γ₀)·M_h_total_eff and its 4×4 eigendecomposition lazily, exposing eigenvalues +
/// eigenvectors as inspectable children.
///
/// <para>This is one slice of a Q-sweep — the "rod" axis whose union over Q forms the
/// 3D / 4D structure that the Object Manager exports.</para>
///
/// <para>Algebraic content: F86 Statement 1 says the slowest pair coalesces at <c>Q_EP =
/// 2/g_eff</c>, where two of the four eigenvalues collide. This snapshot's eigenvalues
/// vector is the local witness — at <c>Q ≈ Q_EP</c> the imaginary parts depart from zero;
/// before/after the algebra is purely real.</para>
/// </summary>
public sealed class LEffSpectralAtQ : IInspectable
{
    public double Q { get; }
    public double J { get; }
    public BlockMatrixIn4Mode LEffView { get; }

    private readonly Lazy<Evd<Complex>> _evd;
    private readonly Lazy<double> _maxImaginaryPart;

    public ComplexVector Eigenvalues => _evd.Value.EigenValues;
    public ComplexMatrix Eigenvectors => _evd.Value.EigenVectors;

    /// <summary>Largest |Im(λ)| among the four eigenvalues. Spikes when Q crosses the EP and
    /// the slowest pair acquires an imaginary part — an empirical EP locator.</summary>
    public double MaxImaginaryPart => _maxImaginaryPart.Value;

    public LEffSpectralAtQ(double q, double gammaZero, ComplexMatrix lEff)
    {
        if (lEff.RowCount != 4 || lEff.ColumnCount != 4)
            throw new ArgumentException($"LEffSpectralAtQ expects 4×4 L_eff; got {lEff.RowCount}×{lEff.ColumnCount}.");
        Q = q;
        J = q * gammaZero;
        LEffView = new BlockMatrixIn4Mode($"L_eff(Q={q:F4})", lEff);

        _evd = new Lazy<Evd<Complex>>(() => lEff.Evd());
        _maxImaginaryPart = new Lazy<double>(() =>
        {
            double m = 0;
            foreach (var z in Eigenvalues)
                if (Math.Abs(z.Imaginary) > m) m = Math.Abs(z.Imaginary);
            return m;
        });
    }

    public string DisplayName => $"L_eff @ Q={Q:F4} (J={J:F4})";
    public string Summary => $"‖L‖_F = {LEffView.Frobenius:F4}, max |Im(λ)| = {MaxImaginaryPart:E2}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return LEffView;
            yield return new InspectableNode("eigenvalues",
                summary: $"max |Im| = {MaxImaginaryPart:E2}",
                payload: new InspectablePayload.Vector("eigenvalues", Eigenvalues, FourModeNames.All));
            yield return new InspectableNode("eigenvectors (columns = right eigvecs)",
                summary: $"4×4 in {string.Join(",", FourModeNames.All)} basis",
                payload: new InspectablePayload.MatrixView("eigenvectors", Eigenvectors,
                    FourModeNames.All, FourModeNames.All));
            yield return InspectableNode.RealScalar("max |Im(λ)|", MaxImaginaryPart, "E3");
            yield return InspectableNode.RealScalar("Q", Q, "F4");
            yield return InspectableNode.RealScalar("J = Q·γ₀", J, "F4");
        }
    }

    public InspectablePayload Payload => LEffView.Payload;
}
