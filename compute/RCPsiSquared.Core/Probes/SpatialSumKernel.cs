using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Probes;

/// <summary>F73-style spatial-sum coherence kernel as a quadratic form on the flat block
/// coefficient vector.
///
/// S(t) = Σ_i 2·|(ρ_i(t))_{0,1}|² can be written as ρ(t)† · S_kernel · ρ(t), where the
/// (0, 1) element of the i-th reduced density matrix is a linear functional of the
/// (n, n+1)-block coefficients (the (p, q) pairs differing only at site i, with p_i = 0
/// and q_i = 1). |·|² becomes A_i^† A_i in matrix form, summed over sites.
///
/// The result is Hermitian and positive semi-definite by construction.
/// </summary>
public static class SpatialSumKernel
{
    public static ComplexMatrix Build(CoherenceBlock block)
    {
        BlockBasis basis = block.Basis;
        int N = block.N;
        int n = block.LowerPopcount;
        int Mtot = basis.MTotal;

        // The A_i vector for site i has nonzeros only at (p, q) where p_i = 0, q = p | bit_i,
        // popcount(q) = n+1. Per site there are at most C(N-1, n) such pairs (choose n bits in
        // the remaining N-1 positions). Collect indices once per site, then update only the
        // dense entries in the index×index outer product.
        var sRaw = new Complex[Mtot, Mtot];
        var siteIndices = new List<int>(capacity: 64);

        for (int site = 0; site < N; site++)
        {
            siteIndices.Clear();
            int maskI = 1 << (N - 1 - site);

            foreach (int p in basis.StatesP)
            {
                if (((p >> (N - 1 - site)) & 1) != 0) continue; // need p_site = 0
                int q = p | maskI;
                if (BitOperations.PopCount((uint)q) != n + 1) continue;
                siteIndices.Add(basis.FlatIndex(p, q));
            }

            // s += 2 · A_i^† A_i on the index×index sub-block. A_i is real-valued (entries 0/1),
            // so the contribution per (i, j) ∈ siteIndices² is +2.
            foreach (int i in siteIndices)
                foreach (int j in siteIndices)
                    sRaw[i, j] += 2.0;
        }

        return Matrix<Complex>.Build.DenseOfArray(sRaw);
    }
}
