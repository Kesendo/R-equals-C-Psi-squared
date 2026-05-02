using System.Text.Json;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Probes;

namespace RCPsiSquared.Cli.Commands;

public static class DecomposeCommand
{
    public static int Run(string[] args)
    {
        var p = new ArgParser(args);
        p.RequireNoPositional();
        int N = p.RequireInt("N");
        int n = p.RequireInt("n");
        double gamma = p.RequireDouble("gamma");
        string? outPath = p.OptionalString("out");

        var block = new CoherenceBlock(N, n, gamma);
        if (block.C < 2)
            throw new InvalidOperationException(
                $"decompose requires chromaticity ≥ 2 (HD=1 and HD=3 channels); got c={block.C} for N={N}, n={n}.");

        var svd = InterChannelSvd.Build(block, hd1: 1, hd2: 3);
        var fourMode = FourModeBasis.Build(block);

        var probe = DickeBlockProbe.Build(block);
        var probe4 = fourMode.Project(probe);

        var result = new
        {
            input = new { N, lower_popcount = n, gamma_0 = gamma },
            chromaticity = block.C,
            block_dim = block.Basis.MTotal,
            sigma_0 = svd.Sigma0,
            top_singular_values = svd.SingularValues.Take(5).ToArray(),
            probe_orthogonal_to_ep_partners = new
            {
                overlap_with_u0 = probe.ConjugateDotProduct(svd.U0InFullBlock).Magnitude,
                overlap_with_v0 = probe.ConjugateDotProduct(svd.V0InFullBlock).Magnitude,
            },
            four_mode_basis = new
            {
                orthonormality_residual = fourMode.OffOrthonormalityResidual,
                probe_components = new
                {
                    c1 = probe4[0].Magnitude,
                    c3 = probe4[1].Magnitude,
                    u0 = probe4[2].Magnitude,
                    v0 = probe4[3].Magnitude,
                },
            },
            ep_naive_q = 2.0 / svd.Sigma0,
        };

        var json = JsonSerializer.Serialize(result, new JsonSerializerOptions { WriteIndented = true });
        if (outPath is null) Console.WriteLine(json);
        else File.WriteAllText(outPath, json);
        return 0;
    }
}
