using System.Globalization;
using System.Text;
using System.Text.Json;
using RCPsiSquared.Core.F89PathK;

namespace RCPsiSquared.Cli.Commands;

/// <summary>Deterministic exact N=6 Route B coefficient export.</summary>
public static class RouteBN6ResidualCommand
{
    public static int Run(string[] args)
    {
        if (args.Length != 2 || args[0] != "--out" || string.IsNullOrWhiteSpace(args[1])
            || args[1].StartsWith("--", StringComparison.Ordinal))
            throw new ArgumentException("usage: route-b-n6-residual --out <path>");

        var even = FoldResultantCertificate.ExportRouteBN6ExactPencil(false);
        var odd = FoldResultantCertificate.ExportRouteBN6ExactPencil(true);
        var payload = new
        {
            schemaVersion = 1,
            n = 6,
            lambdaConvention = "Lambda=2*lambda",
            parameterConvention = "t=i*qCSharp",
            coefficientOrder = "lambda-lowest-first,t-lowest-first",
            rEven = ParityPayload(even),
            rOdd = ParityPayload(odd)
        };
        File.WriteAllText(args[1], JsonSerializer.Serialize(payload) + "\n", new UTF8Encoding(false));
        return 0;
    }

    private static object ParityPayload(RouteBN6ExactPencil pencil) => new
    {
        rOdd = pencil.ROdd,
        sectorDimension = pencil.SectorDimension,
        residualInT = pencil.ResidualInT.Select(row => row.Select(c => c.ToString(CultureInfo.InvariantCulture)).ToArray()).ToArray(),
        atFactorInT = pencil.AtFactorInT.Select(row => row.Select(c => c.ToString(CultureInfo.InvariantCulture)).ToArray()).ToArray()
    };
}
