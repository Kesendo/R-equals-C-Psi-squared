using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.DZero;
using RCPsiSquared.Diagnostics.Foundation;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Lock for the F88-Lens hardware result on IBM Marrakesh framework_snapshots
/// (2026-04-26, job d7mt7jbaq2pc73a24220). The snapshot D F87-trichotomy categories
/// (truly / soft / hard, |+−+⟩ at t = 0.8 with bilinear Hamiltonians, 9-Pauli
/// tomography on (q0, q2)) reconstruct to 2-qubit reduced ρ that, when run through
/// <see cref="MemoryAxisRho"/>, give Π²-odd-fraction-within-memory values
/// differentiating the trichotomy by ~25× between truly and soft. Source data
/// transcribed from
/// <c>D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography\results\framework_snapshots_ibm_marrakesh_20260426_105948.json</c>;
/// Python counterpart: <c>simulations/_f88_lens_ibm_framework_snapshots.py</c>.
/// </summary>
public class IbmFrameworkSnapshotsTrichotomyTests
{
    [Theory]
    [InlineData("truly", 0.0297, 0.5762, 0.4238)]
    [InlineData("soft",  0.7444, 0.5153, 0.4847)]
    [InlineData("hard",  0.2763, 0.3965, 0.6035)]
    public void MarrakeshSnapshotD_F87Trichotomy_Pi2OddInMemory_LocksHardwareValue(
        string category, double expectedPi2Odd, double expectedStatic, double expectedMemory)
    {
        var rho = ReconstructTwoQubitRho(SnapshotDExpectations(category));

        // 2-qubit subsystem reduced from a 3-qubit chain: kernel of L is span{P_0, P_1, P_2}
        // for Heisenberg + Z-dephasing at N = 2 (verified equivalent to popcount-sector projectors).
        var chain = new ChainSystem(N: 2, J: 1.0, GammaZero: 0.05);
        var sm = StationaryModes.Compute(chain);
        var reading = MemoryAxisRho.Decompose(rho, chain, sm);

        Assert.Equal(expectedPi2Odd, reading.Pi2OddFractionWithinMemory, 4);
        Assert.Equal(expectedStatic, reading.StaticFraction, 4);
        Assert.Equal(expectedMemory, reading.MemoryFraction, 4);
    }

    [Fact]
    public void MarrakeshSnapshotD_TrichotomySeparation_TrulyVsSoft_AtLeast20x()
    {
        double trulyOdd = ComputePi2OddInMemory(SnapshotDExpectations("truly"));
        double softOdd = ComputePi2OddInMemory(SnapshotDExpectations("soft"));
        Assert.True(softOdd / trulyOdd > 20.0,
            $"truly→soft separation: soft={softOdd:F4}, truly={trulyOdd:F4}, ratio={softOdd / trulyOdd:F2}");
    }

    private static double ComputePi2OddInMemory(Dictionary<string, double> expectations)
    {
        var rho = ReconstructTwoQubitRho(expectations);
        var chain = new ChainSystem(N: 2, J: 1.0, GammaZero: 0.05);
        var sm = StationaryModes.Compute(chain);
        return MemoryAxisRho.Decompose(rho, chain, sm).Pi2OddFractionWithinMemory;
    }

    private static ComplexMatrix ReconstructTwoQubitRho(Dictionary<string, double> expectations)
    {
        // ρ = (1/4) · Σ_{α, β} ⟨σ_α ⊗ σ_β⟩ · σ_α ⊗ σ_β  (standard inverse Pauli decomposition)
        var rho = ComplexMatrix.Build.Dense(4, 4);
        foreach (var (key, expValue) in expectations)
        {
            string[] parts = key.Split(',');
            var letters = new[] { PauliLetterExtensions.FromSymbol(parts[0][0]), PauliLetterExtensions.FromSymbol(parts[1][0]) };
            var sigma = PauliString.Build(letters);
            rho += (Complex)(expValue / 4.0) * sigma;
        }
        return rho;
    }

    /// <summary>16 Pauli expectation values from snapshot D (q0, q2 reduced after partial
    /// trace over q1). Values verbatim from
    /// framework_snapshots_ibm_marrakesh_20260426_105948.json,
    /// expectations_per_category[category].</summary>
    private static Dictionary<string, double> SnapshotDExpectations(string category) => category switch
    {
        "truly" => new Dictionary<string, double>
        {
            ["I,I"] =  1.0,
            ["X,X"] =  0.00537109375,
            ["X,Y"] = -0.001953125,
            ["X,Z"] =  0.0732421875,
            ["X,I"] =  0.2119140625,
            ["Y,X"] =  0.07470703125,
            ["Y,Y"] =  0.5166015625,
            ["Y,Z"] =  0.62548828125,
            ["Y,I"] =  0.07666015625,
            ["Z,X"] =  0.03662109375,
            ["Z,Y"] =  0.10498046875,
            ["Z,Z"] =  0.1630859375,
            ["Z,I"] =  0.09814453125,
            ["I,X"] =  0.14990234375,
            ["I,Y"] = -0.01513671875,
            ["I,Z"] =  0.00732421875,
        },
        "soft" => new Dictionary<string, double>
        {
            ["I,I"] =  1.0,
            ["X,X"] =  0.00634765625,
            ["X,Y"] = -0.01025390625,
            ["X,Z"] = -0.75146484375,
            ["X,I"] =  0.1435546875,
            ["Y,X"] =  0.0244140625,
            ["Y,Y"] =  0.48779296875,
            ["Y,Z"] = -0.01513671875,
            ["Y,I"] =  0.0107421875,
            ["Z,X"] = -0.49365234375,
            ["Z,Y"] =  0.02880859375,
            ["Z,Z"] =  0.2587890625,
            ["Z,I"] = -0.2333984375,
            ["I,X"] =  0.1376953125,
            ["I,Y"] = -0.03564453125,
            ["I,Z"] = -0.1943359375,
        },
        "hard" => new Dictionary<string, double>
        {
            ["I,I"] =  1.0,
            ["X,X"] =  0.2734375,
            ["X,Y"] =  0.40185546875,
            ["X,Z"] =  0.2255859375,
            ["X,I"] =  0.97265625,
            ["Y,X"] =  0.015625,
            ["Y,Y"] =  0.033203125,
            ["Y,Z"] =  0.01953125,
            ["Y,I"] =  0.02099609375,
            ["Z,X"] = -0.01416015625,
            ["Z,Y"] =  0.015625,
            ["Z,Z"] =  0.01806640625,
            ["Z,I"] =  0.03173828125,
            ["I,X"] =  0.33349609375,
            ["I,Y"] =  0.45166015625,
            ["I,Z"] =  0.2080078125,
        },
        _ => throw new ArgumentException($"unknown trichotomy category: {category}"),
    };
}
