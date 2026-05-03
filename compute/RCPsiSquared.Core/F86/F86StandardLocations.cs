namespace RCPsiSquared.Core.F86;

/// <summary>Single source of truth for the (c, N) sample envelope used by every F86
/// witness/table type. The Python step_f/g sweeps that established the universal-shape
/// Tier-1-candidate ran exactly these points; keeping them in one place means changing the
/// envelope (e.g. adding c=5 once compute supports it) updates all dependent witnesses.
/// </summary>
public static class F86StandardLocations
{
    /// <summary>Full envelope: c=2..4, N=5..8 (with c=4 only at N=7,8 for compute reasons).</summary>
    public static IReadOnlyList<(int C, int N)> Full { get; } = new[]
    {
        (2, 5), (2, 6), (2, 7), (2, 8),
        (3, 5), (3, 6), (3, 7), (3, 8),
        (4, 7), (4, 8),
    };

    /// <summary>Endpoint subset — drops c=4 because the original Endpoint scan at c=4 had
    /// numerical issues at the original sweep grid. Same envelope is fine when scans are
    /// run from C# directly.</summary>
    public static IReadOnlyList<(int C, int N)> EndpointDefault { get; } = new[]
    {
        (2, 5), (2, 6), (2, 7), (2, 8),
        (3, 5), (3, 6), (3, 7), (3, 8),
    };
}
