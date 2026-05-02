using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.DZero;

/// <summary>Kernel of L — the d=0 substrate: eigenvectors of the Liouvillian with |λ| &lt; tol.
///
/// For uniform XY/Heisenberg + Z-dephasing this matches **F4**'s stationary-mode count:
/// N+1 sector projectors P_n on {0, 1, …, N}-excitation subspaces, all living in the
/// {I, Z}^N Pauli sublattice (n_xy = 0). See docs/ANALYTICAL_FORMULAS.md F4 entry.
///
/// This is a verkettung of <see cref="ChainSystem.BuildLiouvillian"/> + MathNet eigendecomp +
/// Pauli-basis projection — no new primitive, just a chained reading.
/// </summary>
public sealed record StationaryModesResult(
    Complex[] Eigenvalues,
    ComplexMatrix RightEigenvectors,
    ComplexMatrix LeftCovectors,
    ComplexMatrix PauliDecomposition,
    int KernelDimension);

public static class StationaryModes
{
    public static StationaryModesResult Compute(ChainSystem chain, ComplexMatrix? Loverride = null, double tolerance = 1e-9)
    {
        var L = Loverride ?? chain.BuildLiouvillian();
        var evd = L.Evd();
        var evals = evd.EigenValues;
        var R = evd.EigenVectors;
        var Rinv = R.Inverse();

        var kernelIdx = new List<int>();
        for (int i = 0; i < evals.Count; i++)
            if (evals[i].Magnitude < tolerance) kernelIdx.Add(i);

        var kEvals = kernelIdx.Select(i => evals[i]).ToArray();
        var kRight = Matrix<Complex>.Build.Dense(R.RowCount, kernelIdx.Count);
        var kLeft = Matrix<Complex>.Build.Dense(kernelIdx.Count, Rinv.ColumnCount);
        for (int i = 0; i < kernelIdx.Count; i++)
        {
            kRight.SetColumn(i, R.Column(kernelIdx[i]));
            kLeft.SetRow(i, Rinv.Row(kernelIdx[i]));
        }

        var transform = PauliBasis.VecToPauliBasisTransform(chain.N);
        var pauliDecomp = (transform.ConjugateTranspose() * kRight) / Math.Pow(2, chain.N);

        return new StationaryModesResult(kEvals, kRight, kLeft, pauliDecomp, kernelIdx.Count);
    }
}
