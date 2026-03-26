using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;

namespace RCPsiSquared.Propagate;

/// <summary>
/// Quantum information metrics on raw Complex[] density matrices.
/// These overloads read from column-major Complex[] arrays without creating
/// MathNet Matrix objects for the full state, enabling N >= 14.
/// Reduced density matrices (1-2 qubits) are small enough for MathNet.
/// </summary>
public static class DensityMatrixToolsRaw
{
    /// <summary>
    /// Partial trace from raw column-major Complex[] to small MathNet Matrix.
    /// Only practical when keep.Length is small (1-2 qubits).
    /// </summary>
    public static Matrix<Complex> PartialTrace(Complex[] rho, int d, int nQubits, int[] keep)
    {
        int nKeep = keep.Length;
        int dKeep = 1 << nKeep;
        var traced = Enumerable.Range(0, nQubits).Except(keep).ToArray();
        int nTraced = traced.Length;

        // Precompute index mapping
        var keptIdx = new int[d];
        var tracedBits = new int[d];
        for (int i = 0; i < d; i++)
        {
            int ki = 0;
            for (int m = 0; m < nKeep; m++)
                ki |= ((i >> (nQubits - 1 - keep[m])) & 1) << (nKeep - 1 - m);
            keptIdx[i] = ki;

            int tb = 0;
            for (int m = 0; m < nTraced; m++)
                tb |= ((i >> (nQubits - 1 - traced[m])) & 1) << (nTraced - 1 - m);
            tracedBits[i] = tb;
        }

        var result = DenseMatrix.Create(dKeep, dKeep, Complex.Zero);

        for (int i = 0; i < d; i++)
        {
            int ki = keptIdx[i];
            int tbi = tracedBits[i];
            for (int j = 0; j < d; j++)
            {
                if (tracedBits[j] != tbi) continue;
                result[ki, keptIdx[j]] += rho[j * d + i]; // column-major
            }
        }
        return result;
    }

    /// <summary>
    /// Mutual information I(A:B) from raw Complex[] density matrix.
    /// </summary>
    public static double MutualInformation(Complex[] rho, int d, int nQubits,
        int[] keepA, int[] keepB)
    {
        var keepAB = keepA.Union(keepB).OrderBy(x => x).ToArray();

        // Trace full state to AB subsystem (small: 4x4 for single-qubit pairs)
        var rhoAB = PartialTrace(rho, d, nQubits, keepAB);

        int nAB = keepAB.Length;
        var localA = keepA.Select(q => Array.IndexOf(keepAB, q)).ToArray();
        var localB = keepB.Select(q => Array.IndexOf(keepAB, q)).ToArray();

        var rhoA = DensityMatrixTools.PartialTrace(rhoAB, nAB, localA);
        var rhoB = DensityMatrixTools.PartialTrace(rhoAB, nAB, localB);

        return DensityMatrixTools.VonNeumannEntropy(rhoA)
             + DensityMatrixTools.VonNeumannEntropy(rhoB)
             - DensityMatrixTools.VonNeumannEntropy(rhoAB);
    }

    /// <summary>
    /// Purity Tr(rho^2) from raw Complex[] array.
    /// </summary>
    public static double Purity(Complex[] rho, int d2)
    {
        double sum = 0;
        // Use parallel reduction for large arrays
        var partialSums = new double[Environment.ProcessorCount];
        Parallel.For(0, d2, () => 0.0,
            (i, _, local) =>
            {
                var v = rho[i];
                return local + v.Real * v.Real + v.Imaginary * v.Imaginary;
            },
            local => { lock (partialSums) { sum += local; } }
        );
        return sum;
    }

    /// <summary>
    /// Build |psi><psi| as raw column-major Complex[] array.
    /// </summary>
    public static Complex[] PureState(Complex[] psi)
    {
        int d = psi.Length;
        var rho = new Complex[d * d];
        Parallel.For(0, d, j =>
        {
            var pjc = Complex.Conjugate(psi[j]);
            int jd = j * d;
            for (int i = 0; i < d; i++)
                rho[jd + i] = psi[i] * pjc;
        });
        return rho;
    }

    /// <summary>
    /// CΨ from raw array, tracing to qubits 0 and 1.
    /// </summary>
    public static (double CPsi, double C, double Psi) ComputeCPsi01(Complex[] rho, int d, int nQubits)
    {
        var rho01 = PartialTrace(rho, d, nQubits, new[] { 0, 1 });
        return DensityMatrixTools.ComputeCPsi(rho01);
    }
}
