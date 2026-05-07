using RCPsiSquared.Cli.Commands;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Cli;

public static class Program
{
    public static int Main(string[] args)
    {
        bool mklActive = MathNetSetup.EnsureInitialized();
        Console.Error.WriteLine($"# MKL: {(mklActive ? "active" : "managed fallback")}, MaxDegreeOfParallelism = {Environment.ProcessorCount}");

        if (args.Length == 0 || args[0] is "-h" or "--help" or "help")
        {
            PrintUsage();
            return 0;
        }

        string command = args[0];
        string[] rest = args[1..];

        try
        {
            return command switch
            {
                "scan" => ScanCommand.Run(rest),
                "decompose" => DecomposeCommand.Run(rest),
                "ep" => EpCommand.Run(rest),
                "plot" => PlotCommand.Run(rest),
                "inspect" => InspectCommand.Run(rest),
                "query" => QueryCommand.Run(rest),
                "knowledge" => KnowledgeCommand.Run(rest),
                _ => UnknownCommand(command),
            };
        }
        catch (Exception ex)
        {
            Console.Error.WriteLine($"error: {ex.Message}");
            return 1;
        }
    }

    private static int UnknownCommand(string command)
    {
        Console.Error.WriteLine($"unknown command: {command}");
        PrintUsage();
        return 2;
    }

    private static void PrintUsage()
    {
        Console.WriteLine("rcpsi: RCPsiSquared.Core CLI");
        Console.WriteLine();
        Console.WriteLine("usage: rcpsi <command> [options]");
        Console.WriteLine();
        Console.WriteLine("commands:");
        Console.WriteLine("  scan        run a resonance scan and print Q_peak / HWHM / |K|max per bond class");
        Console.WriteLine("              args: --N <int> --n <int> --gamma <double> [--out <path>]");
        Console.WriteLine();
        Console.WriteLine("  decompose   compute SVD of inter-channel coupling and 4-mode basis diagnostics");
        Console.WriteLine("              args: --N <int> --n <int> --gamma <double> [--out <path>]");
        Console.WriteLine();
        Console.WriteLine("  ep          print EP-algebra constants t_peak and Q_EP");
        Console.WriteLine("              args: --gamma <double> [--g-eff <double>]");
        Console.WriteLine();
        Console.WriteLine("  plot        run a scan and emit a PNG plot");
        Console.WriteLine("              args: <kind> --N <int> --n <int> --gamma <double> --out <path>");
        Console.WriteLine("              kind: kcurve | shape");
        Console.WriteLine();
        Console.WriteLine("  inspect     walk an IInspectable tree and print it (Object Manager)");
        Console.WriteLine("              args: --N <int> --n <int> --gamma <double> [--root fourmode|f86] [--max-depth 4]");
        Console.WriteLine("                    [--q-sweep] [--with-measured] [--q-grid-points N]");
        Console.WriteLine("                    [--export-json <path>] [--json-only]");
        Console.WriteLine();
        Console.WriteLine("  query       ask typed questions of the F86 knowledge graph");
        Console.WriteLine("              args: --N <int> --n <int> --gamma <double> --q <query> [extra args]");
        Console.WriteLine("              queries: tier-inventory | tier=<name> | anchors |");
        Console.WriteLine("                       witnesses-at --c <int> --wN <int> | per-block-qpeak --c <int> |");
        Console.WriteLine("                       per-bond-qpeak --c <int> --N <int> --bond <Endpoint|Interior> |");
        Console.WriteLine("                       extract-gamma --j <double> --c <int> | open | retracted |");
        Console.WriteLine("                       compare [--q-grid-points N]");
        Console.WriteLine();
        Console.WriteLine("  knowledge   query the typed-knowledge registry built from Runtime + Orchestration");
        Console.WriteLine("              args: <sub> [args]");
        Console.WriteLine("              subs: tier <T1D|T1C|T2V|T2E|OQ|R> | ancestors <Name> | descendants <Name> | all");
        Console.WriteLine();
        Console.WriteLine("examples:");
        Console.WriteLine("  rcpsi scan --N 5 --n 1 --gamma 0.05");
        Console.WriteLine("  rcpsi decompose --N 7 --n 1 --gamma 0.05 --out decomp.json");
        Console.WriteLine("  rcpsi ep --gamma 0.05 --g-eff 2.8284");
        Console.WriteLine("  rcpsi plot kcurve --N 5 --n 1 --gamma 0.05 --out kcurve.png");
        Console.WriteLine("  rcpsi inspect --N 5 --n 1 --gamma 0.05 --max-depth 3");
    }
}
