namespace RCPsiSquared.Core.Confirmations;

/// <summary>A single hardware-confirmed framework prediction. Each entry records a
/// specific prediction-vs-measurement comparison from an IBM Quantum (or similar) run,
/// the framework primitive that produced the prediction, and pointers to the raw data
/// and the experiment writeup.
/// </summary>
public sealed record Confirmation(
    string Name,
    string Date,
    string Machine,
    string JobId,
    string Observable,
    string PredictedValue,
    string MeasuredValue,
    string HardwareData,
    string ExperimentDoc,
    string FrameworkPrimitive,
    string Description);
