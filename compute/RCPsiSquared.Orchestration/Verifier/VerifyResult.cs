using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Orchestration.Verifier;

public sealed record VerifyResult(
    int Total,
    int Drifted,
    IReadOnlyList<DriftReport> Reports
);
