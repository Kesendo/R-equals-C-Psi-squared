using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>HD-channel-uniform projectors for the (n, n+1) coherence block. Generalises
/// **F73** (c=1 spatial-sum closure at the (0, 1) block) to all chromaticities — see
/// docs/ANALYTICAL_FORMULAS.md F73 entry.
///
/// <para>For each channel k ∈ 0..c−1 the channel-uniform vector |c_k⟩ is the equal-weight
/// superposition over all (p, q) basis pairs with HD(p, q) = 2k+1, unit-normalised.</para>
///
/// <para>Empirical fact (verified across N=3..9, multiple topologies, c ∈ {1..4}): the
/// total-H projection P† · M_H_total · P is purely diagonal in the channel-uniform basis.
/// Off-diagonals are exactly zero. This is the structural building block for the heuristic
/// 2-level / 4-mode reduction in F86 (the EP physics lives in the orthogonal complement of
/// this subspace — see PROOF_F86_QPEAK.md).</para>
/// </summary>
public sealed class HdChannelBasis
{
    public CoherenceBlock Block { get; }

    /// <summary>Mtot × c orthonormal projector. Columns are channel-uniform vectors.</summary>
    public ComplexMatrix P { get; }

    /// <summary>Hamming-distance values, ascending: {1, 3, …, 2c−1}.</summary>
    public IReadOnlyList<int> HammingDistances { get; }

    public int C => HammingDistances.Count;

    private HdChannelBasis(CoherenceBlock block, ComplexMatrix p, IReadOnlyList<int> hds)
    {
        Block = block;
        P = p;
        HammingDistances = hds;
    }

    public static HdChannelBasis Build(CoherenceBlock block)
    {
        BlockBasis basis = block.Basis;
        int c = block.C;
        var hds = Chromaticity.HammingDistances(block.N, block.LowerPopcount);
        var p = Matrix<Complex>.Build.Sparse(basis.MTotal, c);

        for (int k = 0; k < c; k++)
        {
            int targetHd = hds[k];
            var indices = new List<int>();
            foreach (long pState in basis.StatesP)
            {
                foreach (long qState in basis.StatesQ)
                {
                    if (BitOperations.PopCount((ulong)(pState ^ qState)) == targetHd)
                        indices.Add(basis.FlatIndex(pState, qState));
                }
            }

            if (indices.Count == 0) continue;
            double norm = 1.0 / Math.Sqrt(indices.Count);
            foreach (int i in indices)
                p[i, k] = new Complex(norm, 0.0);
        }

        return new HdChannelBasis(block, p, hds);
    }
}
