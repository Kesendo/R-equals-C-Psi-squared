namespace RCPsiSquared.Core.Calibration;

/// <summary>Per-qubit metrics from an IBM calibration CSV; mirrors the columns in
/// <c>data/ibm_calibration_snapshots/ibm_marrakesh_calibrations_*.csv</c>. T1, T2,
/// readout error, gate errors, operationality, and the directed coupling graph
/// (neighbour → CZ-error / RZZ-error) remain the measured inputs; the R* properties
/// are derived proxy readings.</summary>
public sealed record QubitData(
    int Qubit,
    double T1Us,
    double T2Us,
    double ReadoutError,
    double SxError,
    double PauliXError,
    bool Operational,
    IReadOnlyDictionary<int, double> CzNeighbours,
    IReadOnlyDictionary<int, double> RzzNeighbours)
{
    /// <summary>r = T2/(2*T1).</summary>
    public double RParam => QubitRegime.RParam(T1Us, T2Us);

    /// <summary>Binary R* band; request an explicit epsilon through
    /// <see cref="QubitRegime.Classify(double,double,double)"/> when a near
    /// band is required.</summary>
    public Regime RStarBand => QubitRegime.Classify(T1Us, T2Us);

    /// <summary>True exactly when r &lt; R*.</summary>
    public bool IsBelowRStar => QubitRegime.IsBelowRStar(T1Us, T2Us);
}
