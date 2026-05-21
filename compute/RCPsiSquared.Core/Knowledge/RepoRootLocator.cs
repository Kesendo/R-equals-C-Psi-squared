using System.IO;

namespace RCPsiSquared.Core.Knowledge;

/// <summary>Locates the repository root by walking up from the current process directory
/// until a directory containing <c>CLAUDE.md</c> is found. Single source shared by the
/// registry's anchor-file validation (<c>ClaimRegistryBuilder</c>), the anchor and wiring
/// audit tests, and the calibration-CSV loaders.</summary>
public static class RepoRootLocator
{
    // Lazy<T> guarantees thread-safe single-execution and an observably-published result.
    // Replaces a hand-rolled (_cached, _searched) pair that had a race window where
    // _searched was set true before _cached was assigned, causing parallel xunit runs
    // to observe null and skip repo-root resolution.
    private static readonly Lazy<string?> _root = new(SearchUpForRepoRoot, isThreadSafe: true);

    /// <summary>The repository root, or <c>null</c> if no ancestor directory contains
    /// <c>CLAUDE.md</c> (e.g. a test process running in an isolated temp directory).</summary>
    public static string? Find() => _root.Value;

    /// <summary>The repository root, throwing <see cref="InvalidOperationException"/> when it
    /// cannot be located. For callers that have no meaningful fallback.</summary>
    public static string Require() =>
        Find() ?? throw new InvalidOperationException(
            "could not locate the repository root: no ancestor of " +
            $"{AppContext.BaseDirectory} contains CLAUDE.md.");

    private static string? SearchUpForRepoRoot()
    {
        var dir = new DirectoryInfo(AppContext.BaseDirectory);
        while (dir != null)
        {
            if (File.Exists(Path.Combine(dir.FullName, "CLAUDE.md")))
                return dir.FullName;
            dir = dir.Parent;
        }
        return null;
    }
}
