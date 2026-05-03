using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Data-dump explorations of the Polynomial → Half → Mirror trio chain. Each
/// fact prints a small data table via ITestOutputHelper. Run with:
/// <code>
/// dotnet test ... --filter "FullyQualifiedName~MemoryAxisChainExploration" \
///                 --logger:"console;verbosity=detailed"
/// </code>
/// </summary>
public class MemoryAxisChainExploration
{
    private readonly ITestOutputHelper _output;
    public MemoryAxisChainExploration(ITestOutputHelper output) => _output = output;

    [Fact]
    public void DumpPi2PartitionAcrossN()
    {
        _output.WriteLine("=== Π² eigenvalue partition (F88 bilinear apex 1/2) ===");
        _output.WriteLine($"{"N",3} {"4^N",10} {"even (+1)",12} {"odd (-1)",12} {"apex",6}");
        _output.WriteLine(new string('-', 50));

        for (int N = 2; N <= 6; N++)
        {
            var chain = new ChainSystem(N, J: 1.0, GammaZero: 0.05);
            var result = MemoryAxisChain.Compute(chain);
            long total = 1L << (2 * N);
            _output.WriteLine($"{N,3} {total,10} {result.Pi2EvenStringCount,12} {result.Pi2OddStringCount,12} {(result.BilinearApexHolds ? "✓" : "✗"),6}");
        }
    }

    [Fact]
    public void DumpF80SpectrumAcrossPi2OddHamiltonians_AtN3()
    {
        _output.WriteLine("=== F80 imaginary spectrum, N = 3 (J = 1.0), Π²-odd 2-body ===");
        _output.WriteLine("");

        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var configs = new (string Name, PauliPairBondTerm[] Terms)[]
        {
            ("XY",      new[] { Bond(PauliLetter.X, PauliLetter.Y) }),
            ("YX",      new[] { Bond(PauliLetter.Y, PauliLetter.X) }),
            ("XZ",      new[] { Bond(PauliLetter.X, PauliLetter.Z) }),
            ("ZX",      new[] { Bond(PauliLetter.Z, PauliLetter.X) }),
            ("XY+YX",   new[] { Bond(PauliLetter.X, PauliLetter.Y), Bond(PauliLetter.Y, PauliLetter.X) }),
            ("XZ+ZX",   new[] { Bond(PauliLetter.X, PauliLetter.Z), Bond(PauliLetter.Z, PauliLetter.X) }),
        };

        foreach (var (name, terms) in configs)
        {
            var result = MemoryAxisChain.Compute(chain, terms);
            _output.WriteLine($"  H = {name,-8}");
            _output.WriteLine($"    Spec(M) = {FormatSpectrum(result.F80ImaginarySpectrum)}");
            _output.WriteLine($"    mirror-symmetric: {result.F80SpectrumIsMirrorSymmetric}");
            _output.WriteLine($"    sum of multiplicities: {result.F80ImaginarySpectrum?.Values.Sum() ?? 0} (expected 4^3 = 64)");
            _output.WriteLine("");
        }
    }

    [Fact]
    public void DumpF80SpectrumScalingWithN_ForXY()
    {
        _output.WriteLine("=== F80 spectrum scaling with N for H = XY (J = 1.0) ===");
        _output.WriteLine("");

        var terms = new[] { Bond(PauliLetter.X, PauliLetter.Y) };
        for (int N = 2; N <= 5; N++)
        {
            var chain = new ChainSystem(N, J: 1.0, GammaZero: 0.05);
            var result = MemoryAxisChain.Compute(chain, terms);
            long total = 1L << (2 * N);

            _output.WriteLine($"  N = {N}:");
            _output.WriteLine($"    Spec(M) = {FormatSpectrum(result.F80ImaginarySpectrum)}");
            _output.WriteLine($"    mirror-symmetric: {result.F80SpectrumIsMirrorSymmetric}");
            _output.WriteLine($"    distinct clusters: {result.F80ImaginarySpectrum?.Count ?? 0}");
            _output.WriteLine($"    total multiplicity: {result.F80ImaginarySpectrum?.Values.Sum() ?? 0} (expected 4^N = {total})");
            _output.WriteLine("");
        }
    }

    [Fact]
    public void DumpTrulyAndPi2EvenNonTruly_BehaviourAtN3()
    {
        _output.WriteLine("=== Truly + Π²-even non-truly behaviour, N = 3 ===");
        _output.WriteLine("");

        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);

        var trulyConfigs = new (string Name, PauliPairBondTerm[] Terms)[]
        {
            ("XX",         new[] { Bond(PauliLetter.X, PauliLetter.X) }),
            ("YY",         new[] { Bond(PauliLetter.Y, PauliLetter.Y) }),
            ("ZZ",         new[] { Bond(PauliLetter.Z, PauliLetter.Z) }),
            ("XX+YY+ZZ",   new[] { Bond(PauliLetter.X, PauliLetter.X),
                                    Bond(PauliLetter.Y, PauliLetter.Y),
                                    Bond(PauliLetter.Z, PauliLetter.Z) }),
        };
        _output.WriteLine("  truly Hamiltonians (M = 0 expected, F80 returns trivial):");
        foreach (var (name, terms) in trulyConfigs)
        {
            var result = MemoryAxisChain.Compute(chain, terms);
            _output.WriteLine($"    {name,-12}: Spec(M) = {FormatSpectrum(result.F80ImaginarySpectrum)}");
        }
        _output.WriteLine("");

        var pi2EvenConfigs = new (string Name, PauliPairBondTerm[] Terms)[]
        {
            ("YZ+ZY",      new[] { Bond(PauliLetter.Y, PauliLetter.Z),
                                    Bond(PauliLetter.Z, PauliLetter.Y) }),
        };
        _output.WriteLine("  Π²-even non-truly (out of F80 scope, expect rejection):");
        foreach (var (name, terms) in pi2EvenConfigs)
        {
            try
            {
                var result = MemoryAxisChain.Compute(chain, terms);
                _output.WriteLine($"    {name,-12}: SURPRISING — F80 accepted; Spec(M) = {FormatSpectrum(result.F80ImaginarySpectrum)}");
            }
            catch (ArgumentException ex)
            {
                string trimmed = ex.Message.Length > 80 ? ex.Message.Substring(0, 80) + "..." : ex.Message;
                _output.WriteLine($"    {name,-12}: rejected (\"{trimmed}\")");
            }
        }
    }

    private static PauliPairBondTerm Bond(PauliLetter a, PauliLetter b) => new(a, b);

    private static string FormatSpectrum(IReadOnlyDictionary<double, int>? spec)
    {
        if (spec == null) return "null";
        if (spec.Count == 0) return "{}";
        return "{ " + string.Join(", ", spec.OrderBy(kv => kv.Key).Select(kv => $"{kv.Key:G6}→{kv.Value}")) + " }";
    }
}
