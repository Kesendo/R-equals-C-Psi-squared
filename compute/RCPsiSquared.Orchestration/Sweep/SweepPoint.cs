using RCPsiSquared.Core.Lindblad;

namespace RCPsiSquared.Orchestration.Sweep;

public sealed record SweepPoint(
    int N,
    HamiltonianClass HClass,
    bool ChainOnly,
    double PredictedValue,
    double LiveValue,
    bool Matched
);
