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

    /// <summary>
    /// Get the Pauli string label for a given index (quaternary encoding).
    /// Index s encodes type_q = (s >> 2q) &amp; 3 for qubit q (bit q).
    /// </summary>
    public static string PauliLabel(int index, int nQubits)
    {
        var chars = new char[nQubits];
        int rem = index;
        for (int q = 0; q < nQubits; q++)
        {
            chars[nQubits - 1 - q] = PauliLabels[rem & 3][0];
            rem >>= 2;
        }
        return new string(chars);
    }

    /// <summary>
    /// XY-weight of a Pauli string: number of X or Y factors.
    /// </summary>
    public static int XYWeight(int index, int nQubits)
    {
        int w = 0;
        int rem = index;
        for (int q = 0; q < nQubits; q++)
        {
            int t = rem & 3;
            if (t == 1 || t == 2) w++; // X=1, Y=2
            rem >>= 2;
        }
        return w;
    }

    // Single-qubit Pauli matrix elements: pauliElement[type, row, col]
    // type: 0=I, 1=X, 2=Y, 3=Z
    private static readonly Complex[,,] PauliElements = InitPauliElements();

    private static Complex[,,] InitPauliElements()
    {
        var p = new Complex[4, 2, 2];
        p[0, 0, 0] = 1; p[0, 1, 1] = 1;                                          // I
        p[1, 0, 1] = 1; p[1, 1, 0] = 1;                                          // X
        p[2, 0, 1] = new Complex(0, -1); p[2, 1, 0] = new Complex(0, 1);         // Y
        p[3, 0, 0] = 1; p[3, 1, 1] = -1;                                         // Z
        return p;
    }

    /// <summary>
    /// Project a superoperator eigenvector onto the Pauli basis.
    /// Input: vec of length 4^N in the |i⟩⟨j| basis (row-major: vec[i*d + j] = ρ[i,j]).
    /// Output: coefficients c_s = Tr(P_s† · unvec(v)) / d for each Pauli string s.
    /// </summary>
    public static Complex[] ProjectOntoPauliBasis(Complex[] vec, int nQubits)
    {
        int d = 1 << nQubits;
        int d2 = d * d;
        int numStrings = 1;
        for (int q = 0; q < nQubits; q++) numStrings *= 4;
        var coeffs = new Complex[numStrings];
        double norm = 1.0 / d;

        for (int s = 0; s < numStrings; s++)
        {
            // Extract base-4 digits of s (Pauli type per qubit)
            var types = new int[nQubits];
            int rem = s;
            for (int q = 0; q < nQubits; q++)
            {
                types[q] = rem & 3;
                rem >>= 2;
            }

            Complex sum = 0;
            for (int i = 0; i < d; i++)
            {
                for (int j = 0; j < d; j++)
                {
                    // P_s[i,j] = product of single-qubit elements
                    Complex pij = 1;
                    bool zero = false;
                    for (int q = 0; q < nQubits; q++)
                    {
                        int iq = (i >> q) & 1;
                        int jq = (j >> q) & 1;
                        var elem = PauliElements[types[q], iq, jq];
                        if (elem == Complex.Zero) { zero = true; break; }
                        pij *= elem;
                    }
                    if (!zero)
                        sum += Complex.Conjugate(pij) * vec[i * d + j];
                }
            }
            coeffs[s] = sum * norm;
        }

        return coeffs;
    }
}
