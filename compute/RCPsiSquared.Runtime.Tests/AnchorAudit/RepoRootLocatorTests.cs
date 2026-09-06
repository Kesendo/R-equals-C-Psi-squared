using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Runtime.Tests.AnchorAudit;

public class RepoRootLocatorTests
{
    [Fact]
    public void Find_UsesTrackedMarkersAvailableInCleanExport()
    {
        var root = RepoRootLocator.Find();

        Assert.NotNull(root);
        Assert.True(File.Exists(Path.Combine(root, "docs", "EXCLUSIONS.md")));
        Assert.True(File.Exists(Path.Combine(root, "compute", "RCPsiSquared.Core",
            "RCPsiSquared.Core.csproj")));
    }

    [Fact]
    public void FindFrom_DoesNotRequireIgnoredClaudeFile()
    {
        var tempRoot = NewTempRoot();
        try
        {
            Directory.CreateDirectory(Path.Combine(tempRoot, "docs"));
            Directory.CreateDirectory(Path.Combine(tempRoot, "compute", "RCPsiSquared.Core"));
            Directory.CreateDirectory(Path.Combine(tempRoot, "a", "b"));
            File.WriteAllText(Path.Combine(tempRoot, "docs", "EXCLUSIONS.md"), "tracked");
            File.WriteAllText(Path.Combine(tempRoot, "compute", "RCPsiSquared.Core",
                "RCPsiSquared.Core.csproj"), "<Project />");

            Assert.Equal(tempRoot, RepoRootLocator.FindFrom(Path.Combine(tempRoot, "a", "b")));
        }
        finally
        {
            Directory.Delete(tempRoot, recursive: true);
        }
    }

    [Fact]
    public void FindFrom_RetainsLegacyClaudeCompatibility()
    {
        var tempRoot = NewTempRoot();
        try
        {
            Directory.CreateDirectory(Path.Combine(tempRoot, "nested"));
            File.WriteAllText(Path.Combine(tempRoot, "CLAUDE.md"), "legacy");

            Assert.Equal(tempRoot, RepoRootLocator.FindFrom(Path.Combine(tempRoot, "nested")));
        }
        finally
        {
            Directory.Delete(tempRoot, recursive: true);
        }
    }

    [Fact]
    public void FindFrom_PrefersTrackedRootOverNestedLegacyMarker()
    {
        var tempRoot = NewTempRoot();
        try
        {
            Directory.CreateDirectory(Path.Combine(tempRoot, "docs"));
            Directory.CreateDirectory(Path.Combine(tempRoot, "compute", "RCPsiSquared.Core"));
            var nested = Path.Combine(tempRoot, "nested");
            Directory.CreateDirectory(nested);
            File.WriteAllText(Path.Combine(tempRoot, "docs", "EXCLUSIONS.md"), "tracked");
            File.WriteAllText(Path.Combine(tempRoot, "compute", "RCPsiSquared.Core",
                "RCPsiSquared.Core.csproj"), "<Project />");
            File.WriteAllText(Path.Combine(nested, "CLAUDE.md"), "legacy nested marker");

            Assert.Equal(tempRoot, RepoRootLocator.FindFrom(nested));
        }
        finally
        {
            Directory.Delete(tempRoot, recursive: true);
        }
    }

    private static string NewTempRoot() =>
        Path.Combine(Path.GetTempPath(), $"rcpsi-root-locator-{Guid.NewGuid():N}");
}
