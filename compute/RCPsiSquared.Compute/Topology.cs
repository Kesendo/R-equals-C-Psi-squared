using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;

namespace RCPsiSquared.Compute;

public record struct Bond(int QubitA, int QubitB, double J, string[] PauliTypes);

/// <summary>
/// Builds Hamiltonians for different quantum network topologies.
/// </summary>
public static class Topology
{
    private static readonly Dictionary<string, Matrix<Complex>> PauliMap = new()
    {
        ["X"] = PauliOps.X,
        ["Y"] = PauliOps.Y,
        ["Z"] = PauliOps.Z,
    };

    public static Bond[] Star(int nQubits, double[] couplings)
    {
        var bonds = new Bond[nQubits - 1];
        for (int i = 1; i < nQubits; i++)
            bonds[i - 1] = new Bond(0, i, couplings[i - 1], new[] { "X", "Y", "Z" });
        return bonds;
    }

    public static Bond[] Chain(int nQubits, double[] couplings)
    {
        var bonds = new Bond[nQubits - 1];
        for (int i = 0; i < nQubits - 1; i++)
            bonds[i] = new Bond(i, i + 1, couplings[i], new[] { "X", "Y", "Z" });
        return bonds;
    }

    /// <summary>
    /// XY chain in the PTF convention: H = Σ (J/2) (X_i X_{i+1} + Y_i Y_{i+1}).
    /// Implemented by halving the bond coupling and restricting PauliTypes to
    /// {X, Y}, so BuildHamiltonian yields (J/2)(XX+YY) per bond.
    /// Used by the PTF dense eigendecomposition (eq014, task PERSPECTIVAL_TIME_FIELD).
    /// </summary>
    public static Bond[] ChainXY(int nQubits, double[] couplings)
    {
        var bonds = new Bond[nQubits - 1];
        for (int i = 0; i < nQubits - 1; i++)
            bonds[i] = new Bond(i, i + 1, couplings[i] / 2.0, new[] { "X", "Y" });
        return bonds;
    }

    public static Bond[] Ring(int nQubits, double[] couplings)
    {
        var bonds = new Bond[nQubits];
        for (int i = 0; i < nQubits; i++)
            bonds[i] = new Bond(i, (i + 1) % nQubits, couplings[i], new[] { "X", "Y", "Z" });
        return bonds;
    }

    public static Bond[] Complete(int nQubits, double J = 1.0)
    {
        var bonds = new List<Bond>();
        for (int i = 0; i < nQubits; i++)
            for (int j = i + 1; j < nQubits; j++)
                bonds.Add(new Bond(i, j, J, new[] { "X", "Y", "Z" }));
        return bonds.ToArray();
    }

    public static Bond[] BinaryTree(int nQubits, double J = 1.0)
    {
        var bonds = new List<Bond>();
        for (int k = 0; k < nQubits; k++)
        {
            int left = 2 * k + 1, right = 2 * k + 2;
            if (left < nQubits) bonds.Add(new Bond(k, left, J, new[] { "X", "Y", "Z" }));
            if (right < nQubits) bonds.Add(new Bond(k, right, J, new[] { "X", "Y", "Z" }));
        }
        return bonds.ToArray();
    }

    /// <summary>
    /// Build the Hamiltonian from bonds.
    /// </summary>
    public static Matrix<Complex> BuildHamiltonian(int nQubits, Bond[] bonds)
    {
        int d = 1 << nQubits;
        Matrix<Complex> H = DenseMatrix.Create(d, d, Complex.Zero);
        foreach (var bond in bonds)
            foreach (var pType in bond.PauliTypes)
                H += bond.J * PauliOps.At(PauliMap[pType], bond.QubitA, nQubits)
                             * PauliOps.At(PauliMap[pType], bond.QubitB, nQubits);
        return H;
    }
}
