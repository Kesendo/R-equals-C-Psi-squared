namespace RCPsiSquared.Core.OpenArcs;

/// <summary>The lifecycle state of an open arc: started research thread that has not been
/// declared finished. <see cref="Open"/> = still parked, awaiting its next step;
/// <see cref="Retired"/> = closed for a recorded reason (subsumed, refuted, or completed
/// elsewhere) and kept only as honest history.</summary>
public enum OpenArcStatus
{
    Open,
    Retired,
}

/// <summary>A single open arc: a research thread that was started, reached a first exemplar,
/// and then parked, with the concrete next move recorded so a later session can pick it up
/// instead of re-discovering that it was never finished.
///
/// <para>The project's chronic failure mode is declaring victory at the first exemplar and
/// forgetting the rest of the arc. This record is the antidote: each entry carries where it
/// stopped (<see cref="ParkedAt"/>) and the next concrete move (<see cref="NextStep"/>), so
/// the world can display its own incompleteness via
/// <see cref="OpenArcsRegistry"/> and the inspect "arcs" section.</para>
/// </summary>
public sealed record OpenArc(
    string Name,
    string Opened,
    string Origin,
    string ParkedAt,
    string NextStep,
    OpenArcStatus Status,
    string? RetiredReason = null);
