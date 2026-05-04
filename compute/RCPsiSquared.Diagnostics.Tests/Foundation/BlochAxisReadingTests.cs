using System.Numerics;
using RCPsiSquared.Core.States;
using RCPsiSquared.Diagnostics.Foundation;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class BlochAxisReadingTests
{
    [Fact]
    public void MaximallyMixed_HasZeroBlochOnAllQubits()
    {
        const int N = 3;
        int d = 1 << N;
        var rhoMm = ComplexMatrix.Build.DiagonalIdentity(d) / d;

        var result = BlochAxisReading.Compute(rhoMm, N);

        Assert.Equal(N, result.Qubits.Count);
        foreach (var q in result.Qubits)
        {
            Assert.Equal(0.0, q.Rx, 12);
            Assert.Equal(0.0, q.Ry, 12);
            Assert.Equal(0.0, q.Rz, 12);
            Assert.Equal(0.0, q.RMagnitude, 12);
            Assert.Equal(0.0, q.EigenDeviation, 12);
            Assert.Equal('I', q.DominantAxis);
            Assert.Equal(0, q.DominantSign);
        }
    }

    [Fact]
    public void PolarityState_PlusUniform_HasXAxisPlusOne_PerQubit()
    {
        // |+++⟩ has each qubit on +X axis: r_k = (1, 0, 0); deviation = 0.5; axis = 'X'.
        const int N = 3;
        var psi = PolarityState.Uniform(N, sign: +1);
        var rho = DensityMatrix.FromStateVector(psi);

        var result = BlochAxisReading.Compute(rho, N);

        foreach (var q in result.Qubits)
        {
            Assert.Equal(1.0, q.Rx, 10);
            Assert.Equal(0.0, q.Ry, 10);
            Assert.Equal(0.0, q.Rz, 10);
            Assert.Equal(1.0, q.RMagnitude, 10);
            Assert.Equal(0.5, q.EigenDeviation, 10);  // ←  the structural ±0.5 pair
            Assert.Equal('X', q.DominantAxis);
            Assert.Equal(+1, q.DominantSign);
        }
    }

    [Fact]
    public void PolarityState_MinusUniform_HasXAxisMinusOne_PerQubit()
    {
        const int N = 3;
        var psi = PolarityState.Uniform(N, sign: -1);
        var rho = DensityMatrix.FromStateVector(psi);

        var result = BlochAxisReading.Compute(rho, N);

        foreach (var q in result.Qubits)
        {
            Assert.Equal(-1.0, q.Rx, 10);
            Assert.Equal('X', q.DominantAxis);
            Assert.Equal(-1, q.DominantSign);
            Assert.Equal(0.5, q.EigenDeviation, 10);  // magnitude is 0.5; sign in DominantSign
        }
    }

    [Fact]
    public void YBasisProduct_HasYAxis_PerQubit()
    {
        // |+i,+i,+i⟩ has each qubit on +Y axis: r_k = (0, 1, 0).
        const int N = 3;
        var psi = YBasisProduct(N, new[] { +1, +1, +1 });
        var rho = DensityMatrix.FromStateVector(psi);

        var result = BlochAxisReading.Compute(rho, N);

        foreach (var q in result.Qubits)
        {
            Assert.Equal(0.0, q.Rx, 10);
            Assert.Equal(1.0, q.Ry, 10);
            Assert.Equal(0.0, q.Rz, 10);
            Assert.Equal(0.5, q.EigenDeviation, 10);
            Assert.Equal('Y', q.DominantAxis);
            Assert.Equal(+1, q.DominantSign);
        }
    }

    [Fact]
    public void ZBasisVacuum_HasZAxisPlusOne_PerQubit()
    {
        // |000⟩ has each qubit on +Z axis: r_k = (0, 0, 1).
        const int N = 3;
        int d = 1 << N;
        var psi = ComplexVector.Build.Dense(d);
        psi[0] = Complex.One;
        var rho = DensityMatrix.FromStateVector(psi);

        var result = BlochAxisReading.Compute(rho, N);

        foreach (var q in result.Qubits)
        {
            Assert.Equal(0.0, q.Rx, 10);
            Assert.Equal(0.0, q.Ry, 10);
            Assert.Equal(1.0, q.Rz, 10);
            Assert.Equal(0.5, q.EigenDeviation, 10);
            Assert.Equal('Z', q.DominantAxis);
            Assert.Equal(+1, q.DominantSign);
        }
    }

    [Fact]
    public void ZBasis_OneOneOne_HasZAxisMinusOne_PerQubit()
    {
        // |111⟩ has each qubit on −Z axis: r_k = (0, 0, -1).
        const int N = 3;
        int d = 1 << N;
        var psi = ComplexVector.Build.Dense(d);
        psi[d - 1] = Complex.One;  // |111⟩ is the last basis state
        var rho = DensityMatrix.FromStateVector(psi);

        var result = BlochAxisReading.Compute(rho, N);

        foreach (var q in result.Qubits)
        {
            Assert.Equal(-1.0, q.Rz, 10);
            Assert.Equal('Z', q.DominantAxis);
            Assert.Equal(-1, q.DominantSign);
        }
    }

    [Fact]
    public void MixedState_001_HasMixedZAxis_PerQubit()
    {
        // |001⟩: q0 = |0⟩ → +Z; q1 = |0⟩ → +Z; q2 = |1⟩ → −Z.
        const int N = 3;
        int d = 1 << N;
        var psi = ComplexVector.Build.Dense(d);
        psi[1] = Complex.One;  // |001⟩ at index 1 (big-endian)
        var rho = DensityMatrix.FromStateVector(psi);

        var result = BlochAxisReading.Compute(rho, N);

        Assert.Equal(+1.0, result.Qubits[0].Rz, 10);
        Assert.Equal(+1.0, result.Qubits[1].Rz, 10);
        Assert.Equal(-1.0, result.Qubits[2].Rz, 10);
    }

    [Fact]
    public void MixedXYZBasis_RecoversAxisPerQubit()
    {
        // |+,+i,0⟩: q0 = X+, q1 = Y+, q2 = Z+.
        const int N = 3;
        var psi = GeneralBasisProduct(N, axes: new[] { 'X', 'Y', 'Z' }, signs: new[] { +1, +1, +1 });
        var rho = DensityMatrix.FromStateVector(psi);

        var result = BlochAxisReading.Compute(rho, N);

        Assert.Equal('X', result.Qubits[0].DominantAxis);
        Assert.Equal('Y', result.Qubits[1].DominantAxis);
        Assert.Equal('Z', result.Qubits[2].DominantAxis);

        foreach (var q in result.Qubits)
        {
            Assert.Equal(+1, q.DominantSign);
            Assert.Equal(0.5, q.EigenDeviation, 10);
        }
    }

    // === Helpers ===

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
                        if ((signs[k] == +1 && bit == 1) || (signs[k] == -1 && bit == 0))
                            zero = true;
                        break;
                }
                if (zero) break;
            }
            vec[idx] = zero ? Complex.Zero : amp;
        }
        return vec;
    }
}
