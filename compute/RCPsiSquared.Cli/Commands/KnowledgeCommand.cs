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
using RCPsiSquared.Runtime.Spectrum;

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
        // Full registry: foundations + Pi2-axes + F88 closed form + F-formula
        // inheritance claims + open questions. Order does not matter (the builder
        // resolves topologically), but grouping reflects the inheritance structure.
        var defaultChain = new ChainSystem(
            N: 5, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.XY, Topology: TopologyKind.Chain);
        var gEff = 1.74;   // pinned Endpoint g_eff from PolarityInheritanceLink

        return new ClaimRegistryBuilder()
            // Foundations
            .RegisterF1Family(defaultChain)
            .RegisterF71Family(N: defaultChain.N)
            .RegisterPi2Family()
            .RegisterF86Main(gammaZero: defaultChain.GammaZero, gEff: gEff)
            .RegisterF86Extended(gammaZero: defaultChain.GammaZero)
            .RegisterF86Item1Light(N: defaultChain.N, n: 1, gammaZero: defaultChain.GammaZero)
            .RegisterHalfIntegerMirror(N: defaultChain.N)
            // Pi2-axes (Halbierungsleiter, Z₄ memory, operator-space mirror)
            .RegisterPi2DyadicLadder()
            .RegisterPi2I4MemoryLoop()
            .RegisterPi2Involution()
            // F88 closed form (memory side anchors)
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .RegisterF88PopcountPairLens(N: defaultChain.N, np: 1, nq: 2)
            // Operator-space mirror (number-side ↔ operator-side per qubit)
            .RegisterPi2OperatorSpaceMirror()
            // Spectrum foundation
            .RegisterW1Dispersion(N: defaultChain.N, J: defaultChain.J, gammaZero: defaultChain.GammaZero)
            // F86 inheritance + meta-claims
            .RegisterF86PolarityLink()
            .RegisterF86LocalGlobalEpLink()
            .RegisterF86PerF71OrbitObservation()
            .RegisterF86JordanWignerLight(N: defaultChain.N, n: 1, gammaZero: defaultChain.GammaZero)
            // F87 family + canonical witnesses
            // RegisterF87Family must follow RegisterPi2Family: DissipatorAxisSelectsPolarityClaim
            // resolves PolarityLayerOriginClaim via b.Get<>() at registration-factory time.
            .RegisterF87Family()
            .RegisterF87StandardWitnessSet()
            // F-formula Pi2-Foundation inheritance claims
            .RegisterF1Pi2Inheritance()
            .RegisterF49Pi2Inheritance()
            .RegisterF80FactorPi2Inheritance()
            .RegisterF81Pi2Inheritance()
            .RegisterF86QEpPi2Inheritance()
            .RegisterF86TPeakPi2Inheritance()
            .RegisterF87Pi2Inheritance()
            .RegisterQubitNecessityPi2Inheritance()
            .RegisterF1T1AmplitudeDampingPi2Inheritance()
            .RegisterF5DepolarizingErrorPi2Inheritance()
            .RegisterDickeSuperpositionQuarterPi2Inheritance()
            .RegisterF39DetPiPi2Inheritance()
            .RegisterF49bCenteredDissipatorPi2Inheritance()
            .RegisterF38Pi2InvolutionPi2Inheritance()
            .RegisterF49cShadowCrossingPi2Inheritance()
            .RegisterF57DwellTimeQuarterPi2Inheritance()
            .RegisterF66PoleModesPi2Inheritance()
            .RegisterF63LCommutesPi2Pi2Inheritance()
            .RegisterF61BitAParityPi2Inheritance()
            .RegisterF77MmSaturationPi2Inheritance()
            .RegisterF60GhzBornBelowFoldPi2Inheritance()
            .RegisterF83AntiFractionPi2Inheritance()
            .RegisterF71MirrorSymmetryPi2Inheritance()
            .RegisterF75MirrorPairMiPi2Inheritance()
            // Open questions
            .RegisterF1OpenQuestions()
            .RegisterF86OpenQuestions()
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
