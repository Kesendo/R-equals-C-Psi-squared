using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;

namespace RCPsiSquared.Propagate;

public record struct Bond(int QubitA, int QubitB, double J, string[] PauliTypes);

/// <summary>
/// Builds Hamiltonians for different quantum network topologies.
/// Copied from RCPsiSquared.Compute, extended with MediatorBridge.
/// </summary>
public static class Topology
{
    private static readonly Dictionary<string, Matrix<Complex>> PauliMap = new()
    {
        ["X"] = PauliOps.X,
        ["Y"] = PauliOps.Y,
        ["Z"] = PauliOps.Z,
    };

    private static readonly string[] Heisenberg = { "X", "Y", "Z" };

    public static Bond[] Chain(int nQubits, double[] couplings)
    {
        var bonds = new Bond[nQubits - 1];
        for (int i = 0; i < nQubits - 1; i++)
            bonds[i] = new Bond(i, i + 1, couplings[i], Heisenberg);
        return bonds;
    }

    public static Bond[] Star(int nQubits, double[] couplings)
    {
        var bonds = new Bond[nQubits - 1];
        for (int i = 1; i < nQubits; i++)
            bonds[i - 1] = new Bond(0, i, couplings[i - 1], Heisenberg);
        return bonds;
    }

    /// <summary>
    /// Recursive mediator bridge topology.
    /// Level 1: 2 qubits (pair)
    /// Level 2: 5 qubits = 2 pairs + mediator
    /// Level 3: 11 qubits = 2 Level-2 bridges + meta-mediator
    /// N(k) = 3 * 2^(k-1) - 1 for k >= 1
    /// </summary>
    public static Bond[] MediatorBridge(int level, double jInternal = 1.0,
        double jBridge = 1.0, double jMeta = 1.0)
    {
        int n = MediatorBridgeSize(level);
        if (level == 1)
            return new[] { new Bond(0, 1, jInternal, Heisenberg) };

        // Level 2: 0-1, 1-2, 2-3, 3-4 (5 qubits)
        if (level == 2)
            return Chain(5, new[] { jInternal, jBridge, jBridge, jInternal });

        // Level 3+: two sub-bridges connected through meta-mediator
        int subSize = MediatorBridgeSize(level - 1);
        int metaQubit = subSize;  // the meta-mediator sits between
        var bonds = new List<Bond>();

        // Left sub-bridge (qubits 0 .. subSize-1)
        var leftBonds = MediatorBridge(level - 1, jInternal, jBridge, jMeta);
        bonds.AddRange(leftBonds);

        // Connect left edge to meta-mediator
        bonds.Add(new Bond(subSize - 1, metaQubit, jMeta, Heisenberg));

        // Connect meta-mediator to right edge
        bonds.Add(new Bond(metaQubit, metaQubit + 1, jMeta, Heisenberg));

        // Right sub-bridge (qubits metaQubit+1 .. n-1)
        var rightBonds = MediatorBridge(level - 1, jInternal, jBridge, jMeta);
        foreach (var b in rightBonds)
            bonds.Add(new Bond(b.QubitA + metaQubit + 1, b.QubitB + metaQubit + 1,
                b.J, b.PauliTypes));

        return bonds.ToArray();
    }

    public static int MediatorBridgeSize(int level)
    {
        if (level <= 0) return 1;
        return 3 * (1 << (level - 1)) - 1;  // 3 * 2^(k-1) - 1
    }

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
