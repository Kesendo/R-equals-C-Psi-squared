namespace RCPsiSquared.Orchestration.Sweep;

public sealed record SweepPoint(
    IReadOnlyDictionary<string, object> Parameters,
    double PredictedValue,
    double LiveValue,
    bool Matched
);
