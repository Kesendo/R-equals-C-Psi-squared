using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Decomposition.Views;

/// <summary>F86 inspection wrapper for D_eff (the dephasing diagonal projected onto the
/// 4-mode basis). Computes the rate-on-each-mode, the trace, and the off-diagonal residual
/// (which must be ~0 for c=2).
///
/// <para>Algebraic content: the 4-mode basis was constructed precisely so that D is diagonal
/// in it; this wrapper is the OOP form of that fact. The off-diagonal residual exposed here
/// is the same quantity tested by <c>FourModeEffectiveTests.Build_DEff_DiagonalRates_AtC2</c>.
/// </para>
/// </summary>
public sealed class DiagonalRatesIn4Mode : IInspectable
{
    public ComplexMatrix Matrix { get; }
    public Complex RateOnC1 { get; }
    public Complex RateOnC3 { get; }
    public Complex RateOnU0 { get; }
    public Complex RateOnV0 { get; }

    private readonly Lazy<double> _offDiagonalResidual;
    private readonly Lazy<Complex> _trace;

    /// <summary>Maximum |off-diagonal element| of D_eff. ~0 verifies that the 4-mode basis
    /// diagonalises D (algebraic identity at c=2).</summary>
    public double OffDiagonalResidual => _offDiagonalResidual.Value;

    /// <summary>Tr(D_eff) = sum of the four mode rates.</summary>
    public Complex Trace => _trace.Value;

    public DiagonalRatesIn4Mode(ComplexMatrix dEff)
    {
        if (dEff.RowCount != 4 || dEff.ColumnCount != 4)
            throw new ArgumentException($"DiagonalRatesIn4Mode expects 4×4; got {dEff.RowCount}×{dEff.ColumnCount}.");
        Matrix = dEff;
        RateOnC1 = dEff[0, 0];
        RateOnC3 = dEff[1, 1];
        RateOnU0 = dEff[2, 2];
        RateOnV0 = dEff[3, 3];
        _offDiagonalResidual = new Lazy<double>(() => MatrixUtilities.MaxOffDiagonalMagnitude(dEff));
        _trace = new Lazy<Complex>(() => RateOnC1 + RateOnC3 + RateOnU0 + RateOnV0);
    }

    public string DisplayName => "D_eff (diagonal rates)";
    public string Summary =>
        $"rates [{RateOnC1.Real:F3}, {RateOnC3.Real:F3}, {RateOnU0.Real:F3}, {RateOnV0.Real:F3}], off-diag = {OffDiagonalResidual:E2}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode($"rate on {FourModeNames.C1}",
                summary: $"{RateOnC1.Real:F4}",
                payload: new InspectablePayload.Scalar(FourModeNames.C1, RateOnC1));
            yield return new InspectableNode($"rate on {FourModeNames.C3}",
                summary: $"{RateOnC3.Real:F4}",
                payload: new InspectablePayload.Scalar(FourModeNames.C3, RateOnC3));
            yield return new InspectableNode($"rate on {FourModeNames.U0}",
                summary: $"{RateOnU0.Real:F4}",
                payload: new InspectablePayload.Scalar(FourModeNames.U0, RateOnU0));
            yield return new InspectableNode($"rate on {FourModeNames.V0}",
                summary: $"{RateOnV0.Real:F4}",
                payload: new InspectablePayload.Scalar(FourModeNames.V0, RateOnV0));
            yield return InspectableNode.RealScalar("off-diagonal residual", OffDiagonalResidual, "E3");
            yield return new InspectableNode("trace",
                summary: $"{Trace.Real:F4}",
                payload: new InspectablePayload.Scalar("Tr(D_eff)", Trace));
        }
    }

    public InspectablePayload Payload =>
        new InspectablePayload.MatrixView("D_eff", Matrix, FourModeNames.All, FourModeNames.All);
}
