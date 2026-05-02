using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Decomposition;

/// <summary>Orthonormal projector onto the full HD-subspace of a coherence block.
///
/// <para>For a fixed Hamming distance HD ∈ {1, 3, …, 2c−1}, the subspace consists of all
/// (p, q) basis pairs with HD(p, q) = HD. The result is a Mtot × n_HD matrix of canonical
/// basis indicators (one nonzero entry per column).</para>
///
/// <para>Note: this is the FULL HD-subspace (e.g. dim 20 for HD=1 at c=2 N=5), not the
/// 1-dimensional channel-uniform vector inside it. The channel-uniform vector is a specific
/// equal-weight linear combination of these basis indicators.</para>
/// </summary>
public static class HdSubspaceProjector
{
    public static ComplexMatrix Build(CoherenceBlock block, int hammingDistance)
    {
        BlockBasis basis = block.Basis;
        var indices = new List<int>();
        foreach (int p in basis.StatesP)
        {
            foreach (int q in basis.StatesQ)
            {
                if (BitOperations.PopCount((uint)(p ^ q)) == hammingDistance)
                    indices.Add(basis.FlatIndex(p, q));
            }
        }

        var proj = Matrix<Complex>.Build.Sparse(basis.MTotal, indices.Count);
        for (int col = 0; col < indices.Count; col++)
            proj[indices[col], col] = Complex.One;
        return proj;
    }
}
