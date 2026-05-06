namespace RCPsiSquared.Core.Knowledge;

/// <summary>Result of a single drift check. A non-drift report carries
/// <see cref="Magnitude"/> = 0 (or null) and a "OK"-shaped description; a drift report
/// carries the magnitude of the deviation and a description naming both the pinned and
/// live values.</summary>
public sealed record DriftReport(
    Type ClaimType,
    bool IsDrift,
    string Description,
    double? Magnitude
);
