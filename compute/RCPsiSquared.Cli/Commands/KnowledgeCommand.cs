using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Orchestration.Cli;
using RCPsiSquared.Orchestration.Render;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Runtime.F71Family;
using RCPsiSquared.Runtime.F86Main;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Cli.Commands;

public static class KnowledgeCommand
{
    public static int Run(string[] args)
    {
        if (args.Length == 0)
        {
            PrintUsage();
            return 0;
        }

        var registry = BuildRegistry();
        var cli = new KnowledgeCli(registry);

        var sub = args[0];
        var rest = args[1..];

        try
        {
            string? output = sub switch
            {
                "tier" => cli.Render(new KnowledgeQuery.Tier(ParseTier(rest))),
                "ancestors" => cli.Render(new KnowledgeQuery.Ancestors(ResolveType(registry, rest))),
                "descendants" => cli.Render(new KnowledgeQuery.Descendants(ResolveType(registry, rest))),
                "all" => cli.Render(new KnowledgeQuery.All()),
                "render" => new KnowledgeRenderer(registry).Render(),
                _ => null,
            };
            if (output is null)
                return UnknownSub(sub);
            Console.Write(output);
            return 0;
        }
        catch (Exception ex)
        {
            Console.Error.WriteLine($"error: {ex.Message}");
            return 1;
        }
    }

    private static ClaimRegistry BuildRegistry()
    {
        // F1 family + Polarity architecture. Add further Register* calls here as new
        // typed-knowledge families are introduced.
        var defaultChain = new ChainSystem(
            N: 5, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.XY, Topology: TopologyKind.Chain);

        return new ClaimRegistryBuilder()
            .RegisterF1Family(defaultChain)
            .RegisterF71Family(N: defaultChain.N)
            .RegisterF86Main(gammaZero: defaultChain.GammaZero, gEff: 1.0)
            .RegisterF86Extended(gammaZero: defaultChain.GammaZero)
            .RegisterF86Item1Light(N: defaultChain.N, n: 1, gammaZero: defaultChain.GammaZero)
            .RegisterPi2Family()
            .RegisterHalfIntegerMirror(N: defaultChain.N)
            .RegisterF86PolarityLink()
            .RegisterF88PopcountCoherence()
            .RegisterF87Family()
            .Build();
    }

    private static Tier ParseTier(string[] args)
    {
        if (args.Length == 0) throw new ArgumentException("missing tier code; expected one of T1D, T1C, T2V, T2E, OQ, R");
        return args[0].ToUpperInvariant() switch
        {
            "T1D" => Tier.Tier1Derived,
            "T1C" => Tier.Tier1Candidate,
            "T2V" => Tier.Tier2Verified,
            "T2E" => Tier.Tier2Empirical,
            "OQ" => Tier.OpenQuestion,
            "R" => Tier.Retracted,
            _ => throw new ArgumentException($"unknown tier code: {args[0]} (expected T1D, T1C, T2V, T2E, OQ, R)"),
        };
    }

    private static Type ResolveType(ClaimRegistry registry, string[] args)
    {
        if (args.Length == 0) throw new ArgumentException("missing claim type name");
        var name = args[0];
        var match = registry.All().FirstOrDefault(c => c.GetType().Name == name);
        if (match is null)
            throw new ArgumentException($"no registered claim with type name '{name}'");
        return match.GetType();
    }

    private static int UnknownSub(string sub)
    {
        Console.Error.WriteLine($"unknown knowledge subcommand: {sub}");
        PrintUsage();
        return 2;
    }

    private static void PrintUsage()
    {
        Console.WriteLine("knowledge: query the typed-knowledge registry");
        Console.WriteLine("usage: rcpsi knowledge <sub> [args]");
        Console.WriteLine("subcommands:");
        Console.WriteLine("  tier <T1D|T1C|T2V|T2E|OQ|R>     list registered Claims at this Tier");
        Console.WriteLine("  ancestors <ClaimTypeName>        transitive parents of a Claim");
        Console.WriteLine("  descendants <ClaimTypeName>      transitive children of a Claim");
        Console.WriteLine("  all                              list every registered Claim");
        Console.WriteLine("  render                           render the registry as Markdown grouped by Tier");
    }
}
