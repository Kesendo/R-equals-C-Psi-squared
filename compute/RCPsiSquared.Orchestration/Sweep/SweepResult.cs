namespace RCPsiSquared.Orchestration.Sweep;

public sealed record SweepResult(
    IReadOnlyList<SweepPoint> Points,
    int Total,
    int Matched,
    int Mismatched
);
