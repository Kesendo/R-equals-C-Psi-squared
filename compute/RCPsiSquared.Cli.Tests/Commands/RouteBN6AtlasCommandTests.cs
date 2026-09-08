using System.Text;
using System.Text.Json;
using RCPsiSquared.Cli.Commands;

namespace RCPsiSquared.Cli.Tests.Commands;

public sealed class RouteBN6AtlasCommandTests
{
    [Fact]
    [Trait("Category", "ROUTE_B_A2_N6_ATLAS_CLI")]
    public void ExportWritesTheGatedAtlasManifestAsDeterministicJson()
    {
        string path = Path.Combine(Path.GetTempPath(), $"route-b-n6-atlas-{Guid.NewGuid():N}.json");
        try
        {
            Assert.Equal(0, RouteBN6AtlasCommand.Run(["--out", path]));
            byte[] bytes = File.ReadAllBytes(path);
            Assert.False(bytes.Take(3).SequenceEqual(new byte[] { 0xef, 0xbb, 0xbf }));
            Assert.Equal((byte)'\n', bytes[^1]);
            Assert.NotEqual((byte)'\r', bytes[^2]);
            using JsonDocument json = JsonDocument.Parse(Encoding.UTF8.GetString(bytes));
            JsonElement root = json.RootElement;
            Assert.Equal(1, root.GetProperty("schemaVersion").GetInt32());
            Assert.Equal(6, root.GetProperty("n").GetInt32());
            Assert.Equal(266, root.GetProperty("loci").GetArrayLength());
            Assert.Equal(96, root.GetProperty("orbits").GetArrayLength());
            Assert.Equal(266, root.GetProperty("reconciliation").GetProperty("consumedLoci").GetInt32());
            Assert.Equal(File.ReadAllBytes(CommittedArtifactPath()), bytes);
        }
        finally
        {
            if (File.Exists(path)) File.Delete(path);
        }
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_N6_ATLAS_CLI")]
    public void InvalidArgumentsAreRejectedWithoutWriting()
    {
        Assert.Throws<ArgumentException>(() => RouteBN6AtlasCommand.Run([]));
        Assert.Throws<ArgumentException>(() => RouteBN6AtlasCommand.Run(["--out"]));
        Assert.Throws<ArgumentException>(() => RouteBN6AtlasCommand.Run(["--unknown", "x"]));
    }

    private static string CommittedArtifactPath()
    {
        DirectoryInfo? cursor = new(AppContext.BaseDirectory);
        while (cursor is not null)
        {
            string candidate = Path.Combine(cursor.FullName, "simulations", "results",
                "route_b_a2_n6_atlas.json");
            if (File.Exists(candidate)) return candidate;
            cursor = cursor.Parent;
        }

        throw new FileNotFoundException("Cannot locate committed N=6 Route-B atlas manifest.");
    }
}
