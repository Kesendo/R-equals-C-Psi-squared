using System.Text;
using System.Text.Json;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Cli.Commands;

/// <summary>Exports the gated N=6 Route-B A2 discovery-atlas manifest.</summary>
public static class RouteBN6AtlasCommand
{
    private static readonly JsonSerializerOptions JsonOptions = new()
    {
        PropertyNamingPolicy = JsonNamingPolicy.CamelCase,
        WriteIndented = true
    };

    public static int Run(string[] args)
    {
        if (args.Length != 2 || args[0] != "--out" || string.IsNullOrWhiteSpace(args[1])
            || args[1].StartsWith("--", StringComparison.Ordinal))
            throw new ArgumentException("usage: route-b-n6-atlas --out <path>");

        var inventory = RouteBA2N6Inventory.LoadDefault();
        RouteBA2N6AtlasManifest manifest =
            new RouteBA2N6CharacterClassifier(inventory).BuildAtlasManifest();
        File.WriteAllText(args[1], JsonSerializer.Serialize(manifest, JsonOptions) + "\n",
            new UTF8Encoding(false));
        return 0;
    }
}
