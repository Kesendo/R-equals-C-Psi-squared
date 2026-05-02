using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Observables;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.States;

namespace RCPsiSquared.Core.Tests.Observables;

public class PiProtectedObservablesTests
{
    [Fact]
    public void Compute_OnZDephasedXYChain_ProducesProtectedAndActiveLists()
    {
        // For a Z-dephased XY chain on |+−+⟩, the framework predicts a non-empty Π-protected
        // skeleton (4^N − 1 candidate Paulis split between protected and active by L's
        // eigenstructure).
        int N = 3;
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var gamma = Enumerable.Repeat(0.1, N).ToArray();
        var psi = PolarityState.Build(N, signs: new[] { +1, -1, +1 });
        var rho0 = DensityMatrix.FromStateVector(psi);

        var result = PiProtectedObservables.Compute(H, gamma, rho0, N);

        // 4^N − 1 = 63 non-identity Paulis total
        Assert.Equal(63, result.Protected.Count + result.Active.Count);
        // At least some Paulis are Π-protected (algebraic skeleton is non-empty)
        Assert.NotEmpty(result.Protected);
        Assert.NotEmpty(result.Active);
        // All protected entries have small max-cluster contribution
        Assert.All(result.Protected, p => Assert.True(p.MaxClusterContribution < 1e-9));
    }

    [Fact]
    public void Compute_PreservesPauliLabels()
    {
        int N = 2;
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var gamma = Enumerable.Repeat(0.1, N).ToArray();
        var psi = PolarityState.Build(N, signs: new[] { +1, -1 });
        var rho0 = DensityMatrix.FromStateVector(psi);

        var result = PiProtectedObservables.Compute(H, gamma, rho0, N);

        // Expect Pauli labels in standard "IXYZ"-format, length N=2
        Assert.All(result.Protected, p => Assert.Equal(2, p.PauliLabel.Length));
        Assert.All(result.Active, p => Assert.Equal(2, p.PauliLabel.Length));
    }
}
