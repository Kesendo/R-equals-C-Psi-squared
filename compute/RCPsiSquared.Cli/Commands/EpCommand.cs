using System.Text.Json;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Cli.Commands;

public static class EpCommand
{
    public static int Run(string[] args)
    {
        var p = new ArgParser(args);
        p.RequireNoPositional();
        double gamma = p.RequireDouble("gamma");
        double? gEff = p.OptionalDouble("g-eff");

        var result = new Dictionary<string, object>
        {
            ["gamma_0"] = gamma,
            ["t_peak"] = EpAlgebra.TPeak(gamma),
        };
        if (gEff.HasValue)
        {
            double qEp = EpAlgebra.QEp(gEff.Value);
            // At Q = Q_EP the EP coupling J·g_eff = 2γ₀ holds by definition (J = γ₀·Q_EP).
            double jAtEp = 2.0 * gamma;
            var (lp, lm) = EpAlgebra.SlowestPairEigenvalues(gamma, jAtEp, gEff.Value);
            result["g_eff"] = gEff.Value;
            result["q_ep"] = qEp;
            result["eigenvalues_at_ep"] = new[] { lp, lm };
        }

        Console.WriteLine(JsonSerializer.Serialize(result, new JsonSerializerOptions { WriteIndented = true }));
        return 0;
    }
}
