namespace RCPsiSquared.Core.Confirmations;

/// <summary>A single hardware-confirmed framework prediction. Each entry records a
/// specific prediction-vs-measurement comparison from an IBM Quantum (or similar) run,
/// the framework primitive that produced the prediction, and pointers to the raw data
/// and the experiment writeup.
///
/// <para><see cref="QubitPath"/> is the physical-qubit list the run targeted
/// (in chain order). Optional and backfilled where the description / hardware-data
/// pointer makes the path unambiguous; entries without a documented path keep it
/// null. Used by <see cref="ConfirmationsRegistry.ByPath"/> /
/// <see cref="ConfirmationsRegistry.ByMachineAndPath"/> for "what has been
/// confirmed on this path?" lookups before a new submission.</para>
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
    string Description,
    IReadOnlyList<int>? QubitPath = null);
