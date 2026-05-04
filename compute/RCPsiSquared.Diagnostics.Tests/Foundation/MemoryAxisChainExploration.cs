using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.States;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit.Abstractions;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

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

    [Fact]
    public void DumpRhoMemoryReadingForCanonicalStates_AtN3()
    {
        _output.WriteLine("=== MemoryAxisRho readings at N = 3 (J = 1.0, γ₀ = 0.05) ===");
        _output.WriteLine($"{"state",-20} {"static",10} {"memory",10} {"Π²-odd/mem",12}");
        _output.WriteLine(new string('-', 56));

        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        int d = 1 << chain.N;

        // Z-basis static states (kernel of L)
        var rhoMm = ComplexMatrix.Build.DiagonalIdentity(d) / d;
        Dump("ρ_mm = I/d", MemoryAxisRho.Decompose(rhoMm, chain));

        var psi000 = ComplexVector.Build.Dense(d);
        psi000[0] = Complex.One;
        Dump("|000⟩ (Z, P_0)", MemoryAxisRho.Decompose(DensityMatrix.FromStateVector(psi000), chain));

        var psi001 = ComplexVector.Build.Dense(d);
        psi001[1] = Complex.One;
        Dump("|001⟩ (Z, P_1)", MemoryAxisRho.Decompose(DensityMatrix.FromStateVector(psi001), chain));

        // X-basis polarity states (Π²-even content only, since {I, X} are both Π²-even)
        Dump("|+++⟩ (X)", MemoryAxisRho.Decompose(
            DensityMatrix.FromStateVector(PolarityState.Uniform(N: 3, sign: +1)), chain));
        Dump("|+−+⟩ (X)", MemoryAxisRho.Decompose(
            DensityMatrix.FromStateVector(PolarityState.Build(N: 3, signs: new[] { +1, -1, +1 })), chain));
        Dump("|+−−⟩ (X)", MemoryAxisRho.Decompose(
            DensityMatrix.FromStateVector(PolarityState.Build(N: 3, signs: new[] { +1, -1, -1 })), chain));

        // Y-basis polarity states ({I, Y}, where Y is Π²-odd; expect non-zero Π²-odd memory)
        Dump("|+i,+i,+i⟩ (Y)", MemoryAxisRho.Decompose(
            DensityMatrix.FromStateVector(YBasisProduct(N: 3, signs: new[] { +1, +1, +1 })), chain));
        Dump("|+i,−i,+i⟩ (Y)", MemoryAxisRho.Decompose(
            DensityMatrix.FromStateVector(YBasisProduct(N: 3, signs: new[] { +1, -1, +1 })), chain));
        Dump("|−i,+i,−i⟩ (Y)", MemoryAxisRho.Decompose(
            DensityMatrix.FromStateVector(YBasisProduct(N: 3, signs: new[] { -1, +1, -1 })), chain));

        // Mixed-basis (X on edges, Y on middle): partial Π²-odd content
        Dump("|+,+i,+⟩ (X-Y-X)", MemoryAxisRho.Decompose(
            DensityMatrix.FromStateVector(MixedXYProduct(N: 3, axes: new[] { 'X', 'Y', 'X' }, signs: new[] { +1, +1, +1 })), chain));
        Dump("|+i,+,+i⟩ (Y-X-Y)", MemoryAxisRho.Decompose(
            DensityMatrix.FromStateVector(MixedXYProduct(N: 3, axes: new[] { 'Y', 'X', 'Y' }, signs: new[] { +1, +1, +1 })), chain));

        void Dump(string name, MemoryAxisRhoResult r)
        {
            _output.WriteLine($"{name,-20} {r.StaticFraction,10:F4} {r.MemoryFraction,10:F4} {r.Pi2OddFractionWithinMemory,12:F4}");
        }
    }

    [Fact]
    public void DumpYZMixedBasisReadings_AtN3()
    {
        _output.WriteLine("=== Y-Z mixed bases (both Π²-odd letter sources) at N = 3 ===");
        _output.WriteLine($"{"state",-22} {"static",10} {"memory",10} {"Π²-odd/mem",12}");
        _output.WriteLine(new string('-', 58));

        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);

        // Y-only and Z-only baselines for comparison
        Dump("|+i,+i,+i⟩ (Y-Y-Y)", new[] { 'Y', 'Y', 'Y' }, new[] { +1, +1, +1 }, chain);
        Dump("|0,0,1⟩    (Z-Z-Z)", new[] { 'Z', 'Z', 'Z' }, new[] { +1, +1, -1 }, chain);

        // Y-Z mix: Y in middle, Z basis on edges (different popcounts)
        Dump("|0,+i,0⟩ (Z-Y-Z, p=0)", new[] { 'Z', 'Y', 'Z' }, new[] { +1, +1, +1 }, chain);
        Dump("|0,+i,1⟩ (Z-Y-Z, p=1)", new[] { 'Z', 'Y', 'Z' }, new[] { +1, +1, -1 }, chain);
        Dump("|1,+i,1⟩ (Z-Y-Z, p=2)", new[] { 'Z', 'Y', 'Z' }, new[] { -1, +1, -1 }, chain);

        // Z in middle, Y on edges
        Dump("|+i,0,+i⟩ (Y-Z-Y, p=0)", new[] { 'Y', 'Z', 'Y' }, new[] { +1, +1, +1 }, chain);
        Dump("|+i,1,+i⟩ (Y-Z-Y, p=1)", new[] { 'Y', 'Z', 'Y' }, new[] { +1, -1, +1 }, chain);

        // Three-letter mix: X-Y-Z (all axis types in one state)
        Dump("|+,+i,0⟩ (X-Y-Z)", new[] { 'X', 'Y', 'Z' }, new[] { +1, +1, +1 }, chain);
        Dump("|+,+i,1⟩ (X-Y-Z)", new[] { 'X', 'Y', 'Z' }, new[] { +1, +1, -1 }, chain);

        // X-Z (no Y; Z is also Π²-odd, so still expect non-zero)
        Dump("|+,0,+⟩ (X-Z-X, p=0)", new[] { 'X', 'Z', 'X' }, new[] { +1, +1, +1 }, chain);
        Dump("|+,1,+⟩ (X-Z-X, p=1)", new[] { 'X', 'Z', 'X' }, new[] { +1, -1, +1 }, chain);

        void Dump(string name, char[] axes, int[] signs, ChainSystem chainSys)
        {
            var psi = GeneralBasisProduct(N: 3, axes: axes, signs: signs);
            var rho = DensityMatrix.FromStateVector(psi);
            var r = MemoryAxisRho.Decompose(rho, chainSys);
            _output.WriteLine($"{name,-22} {r.StaticFraction,10:F4} {r.MemoryFraction,10:F4} {r.Pi2OddFractionWithinMemory,12:F4}");
        }
    }

    /// <summary>Y-basis polarity tensor product: |+i⟩ = (|0⟩ + i|1⟩)/√2, |−i⟩ = (|0⟩ − i|1⟩)/√2,
    /// per-site signs +1 (=|+i⟩) or −1 (=|−i⟩). Big-endian site convention matches PolarityState.</summary>
    private static ComplexVector YBasisProduct(int N, IReadOnlyList<int> signs)
    {
        int d = 1 << N;
        var vec = ComplexVector.Build.Dense(d);
        double norm = 1.0 / Math.Sqrt(d);

        for (int idx = 0; idx < d; idx++)
        {
            Complex amp = Complex.One;
            for (int k = 0; k < N; k++)
            {
                int bit = (idx >> (N - 1 - k)) & 1;
                if (bit == 1)
                    amp *= signs[k] == +1 ? Complex.ImaginaryOne : -Complex.ImaginaryOne;
            }
            vec[idx] = amp * norm;
        }
        return vec;
    }

    /// <summary>Mixed X/Y-basis tensor product: each site is in either σ_x or σ_y eigenstate
    /// per <paramref name="axes"/> ('X' or 'Y'), with sign +1 (|+⟩ or |+i⟩) or −1 (|−⟩ or |−i⟩).</summary>
    private static ComplexVector MixedXYProduct(int N, IReadOnlyList<char> axes, IReadOnlyList<int> signs)
    {
        int d = 1 << N;
        var vec = ComplexVector.Build.Dense(d);
        double norm = 1.0 / Math.Sqrt(d);

        for (int idx = 0; idx < d; idx++)
        {
            Complex amp = Complex.One;
            for (int k = 0; k < N; k++)
            {
                int bit = (idx >> (N - 1 - k)) & 1;
                if (bit == 1)
                {
                    if (axes[k] == 'X')
                        amp *= signs[k]; // +1 for |+⟩, −1 for |−⟩
                    else // 'Y'
                        amp *= signs[k] == +1 ? Complex.ImaginaryOne : -Complex.ImaginaryOne;
                }
            }
            vec[idx] = amp * norm;
        }
        return vec;
    }

    [Fact]
    public void DumpPerQubitBlochReading_AcrossCanonicalStates_AtN3()
    {
        _output.WriteLine("=== Per-qubit Bloch ±0.5-pair reading at N = 3 ===");
        _output.WriteLine("(EigenDeviation = |r_k|/2; the structural ±-pair around 1/2)");
        _output.WriteLine($"{"state",-22} {"q0 (axis ± dev)",18} {"q1 (axis ± dev)",18} {"q2 (axis ± dev)",18}");
        _output.WriteLine(new string('-', 80));

        const int N = 3;
        int d = 1 << N;

        // Mixed
        var rhoMm = ComplexMatrix.Build.DiagonalIdentity(d) / d;
        Dump("ρ_mm = I/d", rhoMm);

        // Z-basis
        var psi000 = ComplexVector.Build.Dense(d);
        psi000[0] = Complex.One;
        Dump("|000⟩", DensityMatrix.FromStateVector(psi000));

        var psi001 = ComplexVector.Build.Dense(d);
        psi001[1] = Complex.One;
        Dump("|001⟩", DensityMatrix.FromStateVector(psi001));

        var psi111 = ComplexVector.Build.Dense(d);
        psi111[d - 1] = Complex.One;
        Dump("|111⟩", DensityMatrix.FromStateVector(psi111));

        // X-basis
        Dump("|+++⟩", DensityMatrix.FromStateVector(PolarityState.Uniform(N, +1)));
        Dump("|−−−⟩", DensityMatrix.FromStateVector(PolarityState.Uniform(N, -1)));
        Dump("|+−+⟩", DensityMatrix.FromStateVector(
            PolarityState.Build(N, new[] { +1, -1, +1 })));

        // Y-basis
        Dump("|+i,+i,+i⟩", DensityMatrix.FromStateVector(
            GeneralBasisProduct(N, new[] { 'Y', 'Y', 'Y' }, new[] { +1, +1, +1 })));
        Dump("|+i,−i,+i⟩", DensityMatrix.FromStateVector(
            GeneralBasisProduct(N, new[] { 'Y', 'Y', 'Y' }, new[] { +1, -1, +1 })));

        // Mixed
        Dump("|+,+i,0⟩", DensityMatrix.FromStateVector(
            GeneralBasisProduct(N, new[] { 'X', 'Y', 'Z' }, new[] { +1, +1, +1 })));
        Dump("|+,+i,1⟩", DensityMatrix.FromStateVector(
            GeneralBasisProduct(N, new[] { 'X', 'Y', 'Z' }, new[] { +1, +1, -1 })));
        Dump("|0,+i,1⟩ (Z-Y-Z)", DensityMatrix.FromStateVector(
            GeneralBasisProduct(N, new[] { 'Z', 'Y', 'Z' }, new[] { +1, +1, -1 })));

        void Dump(string name, ComplexMatrix rho)
        {
            var r = BlochAxisReading.Compute(rho, N);
            _output.WriteLine($"{name,-22} {Format(r.Qubits[0]),18} {Format(r.Qubits[1]),18} {Format(r.Qubits[2]),18}");
        }

        static string Format(QubitBlochReading q)
        {
            if (q.DominantAxis == 'I') return "I  (mixed)";
            char sign = q.DominantSign > 0 ? '+' : '−';
            return $"{q.DominantAxis} {sign} {q.EigenDeviation:F3}";
        }
    }

    /// <summary>General X/Y/Z-basis tensor product. axes[k] ∈ {'X', 'Y', 'Z'}; signs[k] ∈ {+1, −1}.
    /// For 'X': ±1 = |+⟩/|−⟩; 'Y': ±1 = |+i⟩/|−i⟩; 'Z': ±1 = |0⟩/|1⟩ (computational basis).</summary>
    private static ComplexVector GeneralBasisProduct(int N, IReadOnlyList<char> axes, IReadOnlyList<int> signs)
    {
        int d = 1 << N;
        var vec = ComplexVector.Build.Dense(d);
        double sqrt2 = Math.Sqrt(2.0);

        for (int idx = 0; idx < d; idx++)
        {
            Complex amp = Complex.One;
            bool zero = false;
            for (int k = 0; k < N; k++)
            {
                int bit = (idx >> (N - 1 - k)) & 1;
                switch (axes[k])
                {
                    case 'X':
                        amp *= (bit == 0 ? 1.0 : (double)signs[k]) / sqrt2;
                        break;
                    case 'Y':
                        if (bit == 0) amp *= 1.0 / sqrt2;
                        else amp *= (signs[k] == +1 ? Complex.ImaginaryOne : -Complex.ImaginaryOne) / sqrt2;
                        break;
                    case 'Z':
                        // |0⟩ (sign=+1) keeps bit=0, drops bit=1; |1⟩ (sign=−1) keeps bit=1
                        if ((signs[k] == +1 && bit == 1) || (signs[k] == -1 && bit == 0))
                        {
                            zero = true;
                        }
                        break;
                    default:
                        throw new ArgumentException($"unknown axis '{axes[k]}'; expected X/Y/Z");
                }
                if (zero) break;
            }
            vec[idx] = zero ? Complex.Zero : amp;
        }
        return vec;
    }

    private static PauliPairBondTerm Bond(PauliLetter a, PauliLetter b) => new(a, b);

    private static string FormatSpectrum(IReadOnlyDictionary<double, int>? spec)
    {
        if (spec == null) return "null";
        if (spec.Count == 0) return "{}";
        return "{ " + string.Join(", ", spec.OrderBy(kv => kv.Key).Select(kv => $"{kv.Key:G6}→{kv.Value}")) + " }";
    }
}
