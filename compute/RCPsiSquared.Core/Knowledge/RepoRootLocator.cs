using System.IO;

namespace RCPsiSquared.Core.Knowledge;

/// <summary>Locates the repository root by walking up from the current process directory.
/// The primary marker is the tracked pair <c>docs/EXCLUSIONS.md</c> and
/// <c>compute/RCPsiSquared.Core/RCPsiSquared.Core.csproj</c>; <c>CLAUDE.md</c>
/// remains a compatibility fallback. Single source shared by the
/// registry's anchor-file validation (<c>ClaimRegistryBuilder</c>), the anchor and wiring
/// audit tests, and the calibration-CSV loaders.</summary>
public static class RepoRootLocator
{
    // Lazy<T> guarantees thread-safe single-execution and an observably-published result.
    // Replaces a hand-rolled (_cached, _searched) pair that had a race window where
    // _searched was set true before _cached was assigned, causing parallel xunit runs
    // to observe null and skip repo-root resolution.
    private static readonly Lazy<string?> _root = new(SearchUpForRepoRoot, isThreadSafe: true);

    /// <summary>The repository root, or <c>null</c> if no ancestor contains the tracked
    /// repository markers or the legacy <c>CLAUDE.md</c> marker.</summary>
    public static string? Find() => _root.Value;

    /// <summary>Search from an explicit directory. This makes clean-export behavior
    /// testable without depending on the test host's <see cref="AppContext.BaseDirectory"/>.</summary>
    public static string? FindFrom(string startDirectory)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(startDirectory);
        return SearchUpForRepoRoot(new DirectoryInfo(Path.GetFullPath(startDirectory)));
    }

    /// <summary>The repository root, throwing <see cref="InvalidOperationException"/> when it
    /// cannot be located. For callers that have no meaningful fallback.</summary>
    public static string Require() =>
        Find() ?? throw new InvalidOperationException(
            "could not locate the repository root: no ancestor of " +
            $"{AppContext.BaseDirectory} contains the tracked repository markers " +
            "or legacy CLAUDE.md.");

    private static string? SearchUpForRepoRoot() =>
        SearchUpForRepoRoot(new DirectoryInfo(AppContext.BaseDirectory));

    private static string? SearchUpForRepoRoot(DirectoryInfo start)
    {
        // Prefer tracked markers across the whole ancestor chain. An ignored
        // CLAUDE.md in a nested directory must not shadow the actual root.
        var dir = start;
        while (dir != null)
        {
            if (HasTrackedMarkers(dir.FullName))
                return dir.FullName;
            dir = dir.Parent;
        }

        dir = start;
        while (dir != null)
        {
            if (File.Exists(Path.Combine(dir.FullName, "CLAUDE.md")))
                return dir.FullName;
            dir = dir.Parent;
        }
        return null;
    }

    private static bool HasTrackedMarkers(string directory) =>
        File.Exists(Path.Combine(directory, "docs", "EXCLUSIONS.md")) &&
        File.Exists(Path.Combine(directory, "compute", "RCPsiSquared.Core",
            "RCPsiSquared.Core.csproj"));
}
