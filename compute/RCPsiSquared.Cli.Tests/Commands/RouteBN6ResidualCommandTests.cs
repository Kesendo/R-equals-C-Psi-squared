using System.Text;
using System.Text.Json;
using RCPsiSquared.Cli.Commands;

namespace RCPsiSquared.Cli.Tests.Commands;

[Trait("Category", "ROUTE_B_A2_N6_PRODUCER")]
public sealed class RouteBN6ResidualCommandTests
{
    [Fact]
    public void ExportIsDeterministicExactJsonWithoutBomAndWithOneFinalNewline()
    {
        var first = Path.Combine(Path.GetTempPath(), $"route-b-n6-{Guid.NewGuid():N}.json");
        var second = Path.Combine(Path.GetTempPath(), $"route-b-n6-{Guid.NewGuid():N}.json");
        try
        {
            Assert.Equal(0, RouteBN6ResidualCommand.Run(new[] { "--out", first }));
            Assert.Equal(0, RouteBN6ResidualCommand.Run(new[] { "--out", second }));
            var bytes = File.ReadAllBytes(first);
            Assert.Equal(bytes, File.ReadAllBytes(second));
            Assert.False(bytes.Take(3).SequenceEqual(new byte[] { 0xef, 0xbb, 0xbf }));
            Assert.Equal((byte)'\n', bytes[^1]);
            Assert.NotEqual((byte)'\n', bytes[^2]);
            Assert.NotEqual((byte)'\r', bytes[^2]);
            using var json = JsonDocument.Parse(Encoding.UTF8.GetString(bytes));
            var root = json.RootElement;
            Assert.Equal(1, root.GetProperty("schemaVersion").GetInt32());
            Assert.Equal(6, root.GetProperty("n").GetInt32());
            Assert.Equal("Lambda=2*lambda", root.GetProperty("lambdaConvention").GetString());
            Assert.Equal("t=i*qCSharp", root.GetProperty("parameterConvention").GetString());
            Assert.Equal("lambda-lowest-first,t-lowest-first", root.GetProperty("coefficientOrder").GetString());
            foreach (var (name, rOdd) in new[] { ("rEven", false), ("rOdd", true) })
            {
                var parity = root.GetProperty(name);
                Assert.Equal(rOdd, parity.GetProperty("rOdd").GetBoolean());
                Assert.Equal(45, parity.GetProperty("sectorDimension").GetInt32());
                Assert.Equal(33, parity.GetProperty("residualInT").GetArrayLength());
                Assert.Equal(14, parity.GetProperty("atFactorInT").GetArrayLength());
                foreach (var key in new[] { "residualInT", "atFactorInT" })
                    foreach (var row in parity.GetProperty(key).EnumerateArray())
                        foreach (var coefficient in row.EnumerateArray())
                            Assert.Matches(@"^-?(0|[1-9][0-9]*)$", coefficient.GetString()!);
            }
        }
        finally
        {
            if (File.Exists(first)) File.Delete(first);
            if (File.Exists(second)) File.Delete(second);
        }
    }

    [Fact]
    public void InvalidArgumentsAreRejectedWithoutWriting()
    {
        var path = Path.Combine(Path.GetTempPath(), $"route-b-n6-invalid-{Guid.NewGuid():N}.json");
        foreach (var args in new[]
        {
            Array.Empty<string>(), new[] { "--out" }, new[] { "--out", "" },
            new[] { "--out", " " }, new[] { "--out", "--unknown" },
            new[] { "--out", path, "--out", path }, new[] { "--unknown", path },
            new[] { "--out", path, "--unknown" }, new[] { path }
        })
            Assert.Throws<ArgumentException>(() => RouteBN6ResidualCommand.Run(args));
        Assert.False(File.Exists(path));
    }
}
