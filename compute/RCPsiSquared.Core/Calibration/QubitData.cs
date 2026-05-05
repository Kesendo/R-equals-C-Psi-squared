namespace RCPsiSquared.Core.Calibration;

/// <summary>Per-qubit calibration metrics from an IBM Heron r2 calibration CSV.
/// Mirrors the columns in <c>ClaudeTasks/IBM_R2_calibration_ibm_marrakesh/
/// ibm_marrakesh_calibrations_*.csv</c>: T1, T2, readout error, single-qubit
/// gate errors, and the directed coupling graph encoded as
/// neighbour → CZ-error / RZZ-error dictionaries.
///
/// <para>Derived regime properties (<see cref="RParam"/>, <see cref="Regime"/>,
/// <see cref="IsQuantumSide"/>) bridge the raw calibration to the framework's
/// CΨ = ¼ fold-catastrophe boundary; see <see cref="QubitRegime"/>.</para></summary>
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
    /// <summary>r = T2 / (2·T1); see <see cref="QubitRegime.RParam"/>.</summary>
    public double RParam => QubitRegime.RParam(T1Us, T2Us);

    /// <summary>Binary regime classification (no boundary band); see
    /// <see cref="QubitRegime.Classify(double, double, double)"/> for the ε-band variant.</summary>
    public Regime Regime => QubitRegime.Classify(T1Us, T2Us);

    /// <summary>True iff r &lt; R*: CΨ_min crosses ¼, bridge open in the d=0 ↔ d=2 sense.</summary>
    public bool IsQuantumSide => QubitRegime.IsQuantumSide(T1Us, T2Us);
}
