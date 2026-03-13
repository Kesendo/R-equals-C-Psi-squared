using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;

namespace RCPsiSquared.Compute;

/// <summary>
/// Single-qubit Pauli matrices and tensor product utilities.
/// All matrices are 2x2 complex. Tensor products build N-qubit operators.
/// </summary>
public static class PauliOps
{
    public static readonly Matrix<Complex> I2 = DenseMatrix.OfArray(new Complex[,]
    {
        { 1, 0 },
        { 0, 1 }
    });

    public static readonly Matrix<Complex> X = DenseMatrix.OfArray(new Complex[,]
    {
        { 0, 1 },
        { 1, 0 }
    });

    public static readonly Matrix<Complex> Y = DenseMatrix.OfArray(new Complex[,]
    {
        { 0, new Complex(0, -1) },
        { new Complex(0, 1), 0 }
    });

    public static readonly Matrix<Complex> Z = DenseMatrix.OfArray(new Complex[,]
    {
        { 1, 0 },
        { 0, -1 }
    });

    public static readonly Matrix<Complex> SigmaMinus = DenseMatrix.OfArray(new Complex[,]
    {
        { 0, 0 },
        { 1, 0 }
    });

    private static readonly Matrix<Complex>[] Paulis = { I2, X, Y, Z };
    private static readonly string[] PauliLabels = { "I", "X", "Y", "Z" };

    /// <summary>
    /// Place operator on qubit `target` in an N-qubit system.
    /// Result = I ⊗ ... ⊗ op ⊗ ... ⊗ I
    /// </summary>
    public static Matrix<Complex> At(Matrix<Complex> op, int target, int nQubits)
    {
        Matrix<Complex> result = target == 0 ? op : I2;
        for (int k = 1; k < nQubits; k++)
            result = result.KroneckerProduct(k == target ? op : I2);
        return result;
    }

    /// <summary>
    /// Generate all N-qubit Pauli strings with labels.
    /// Yields (label, matrix) for each of 4^N operators.
    /// </summary>
    public static IEnumerable<(string Label, Matrix<Complex> Op)> PauliBasis(int nQubits)
    {
        int total = (int)Math.Pow(4, nQubits);
        for (int idx = 0; idx < total; idx++)
        {
            string label = "";
            Matrix<Complex> op = DenseMatrix.OfArray(new Complex[,] { { 1 } });
            int rem = idx;
            for (int q = 0; q < nQubits; q++)
            {
                int pIdx = rem % 4;
                rem /= 4;
                label = PauliLabels[pIdx] + label;
                op = Paulis[pIdx].KroneckerProduct(op);
            }
            yield return (label, op);
        }
    }
}
