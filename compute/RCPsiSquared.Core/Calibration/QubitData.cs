namespace RCPsiSquared.Core.Calibration;

/// <summary>Per-qubit calibration metrics from an IBM Heron r2 calibration CSV.
/// Mirrors the columns in <c>ClaudeTasks/IBM_R2_calibration_ibm_marrakesh/
/// ibm_marrakesh_calibrations_*.csv</c>: T1, T2, readout error, single-qubit
/// gate errors, and the directed coupling graph encoded as
/// neighbour → CZ-error / RZZ-error dictionaries.</summary>
public sealed record QubitData(
    int Qubit,
    double T1Us,
    double T2Us,
    double ReadoutError,
    double SxError,
    double PauliXError,
    bool Operational,
    IReadOnlyDictionary<int, double> CzNeighbours,
    IReadOnlyDictionary<int, double> RzzNeighbours
);
