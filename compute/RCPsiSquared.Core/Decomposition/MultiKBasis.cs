using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Decomposition;

/// <summary>F86 Item 2: extension of <see cref="FourModeBasis"/> to chromaticity c ≥ 3.
/// For each adjacent rate-channel pair (HD = 2k−1, HD = 2k+1) with k ∈ {1, …, c−1}, build
/// the quartet (|c_{2k−1}⟩, |c_{2k+1}⟩, |u_0^{(k)}⟩, |v_0^{(k)}⟩) just like the 4-mode
/// basis. Concatenate all c−1 quartets into a single 4(c−1)-dimensional orthonormal basis.
///
/// <para>Slowest-pair-dominance hypothesis: at low Q (around the canonical Interior peak),
/// the k=1 quartet alone reproduces the K-curve. Higher-k quartets contribute only at
/// higher Q. This basis lets us test that hypothesis and recover analytical structure
/// across all c.</para>
/// </summary>
public sealed class MultiKBasis
{
    public CoherenceBlock Block { get; }
    public IReadOnlyList<KQuartet> Quartets { get; }
    public ComplexMatrix BasisMatrix { get; }
    public double OffOrthonormalityResidual { get; }

    public int TotalModes => BasisMatrix.ColumnCount;

    private MultiKBasis(CoherenceBlock block, IReadOnlyList<KQuartet> quartets,
        ComplexMatrix basis, double residual)
    {
        Block = block;
        Quartets = quartets;
        BasisMatrix = basis;
        OffOrthonormalityResidual = residual;
    }

    /// <summary>Build the multi-k basis. Each k ∈ {1, …, c−1} provides a quartet; the
    /// channel-uniform vector |c_{2k+1}⟩ is shared between adjacent quartets (it's both the
    /// "+1" of pair k and the "−1" of pair k+1), so we deduplicate. The c−2 shared
    /// channel-uniforms plus possible non-orthogonal SVD-top vectors across HD subspaces
    /// mean the full candidate set is over-complete; we run Gram-Schmidt with a tolerance
    /// to extract an orthonormal subset of rank ≤ 3c−2.</summary>
    public static MultiKBasis Build(CoherenceBlock block, double orthoTolerance = 1e-10)
    {
        if (block.C < 2)
            throw new ArgumentException($"MultiKBasis requires chromaticity ≥ 2; got c={block.C}.");

        var hdChannel = HdChannelBasis.Build(block);
        int c = block.C;
        int dim = block.Basis.MTotal;

        // Step 1: collect candidate vectors in order — all channel-uniform vectors first,
        // then SVD top vectors per quartet.
        var quartets = new KQuartet[c - 1];
        var candidates = new List<ComplexVector>();

        // All c channel-uniform vectors (|c_1⟩, |c_3⟩, …, |c_{2c−1}⟩)
        var cuVectors = new ComplexVector[c];
        for (int k = 0; k < c; k++)
        {
            int hd = 2 * k + 1;
            int idx = -1;
            for (int i = 0; i < hdChannel.HammingDistances.Count; i++)
                if (hdChannel.HammingDistances[i] == hd) { idx = i; break; }
            if (idx < 0)
                throw new InvalidOperationException($"HD={hd} not present at c={block.C}");
            cuVectors[k] = hdChannel.P.Column(idx);
            candidates.Add(cuVectors[k]);
        }

        // SVD vectors per k
        for (int k = 1; k <= c - 1; k++)
        {
            int hd1 = 2 * k - 1;
            int hd2 = 2 * k + 1;
            var svd = InterChannelSvd.Build(block, hd1, hd2);
            quartets[k - 1] = new KQuartet(k, hd1, hd2, svd.Sigma0);
            candidates.Add(svd.U0InFullBlock);
            candidates.Add(svd.V0InFullBlock);
        }

        // Step 2: Gram-Schmidt orthonormalisation, dropping vectors that fall below tolerance.
        // (No defensive clone needed — operator `-` returns a new vector each iteration, so
        // the original candidate is never mutated.)
        var orthonormal = new List<ComplexVector>();
        foreach (var candidate in candidates)
        {
            var v = candidate;
            foreach (var u in orthonormal)
                v -= u.ConjugateDotProduct(v) * u;
            double norm = v.L2Norm();
            if (norm > orthoTolerance)
                orthonormal.Add(v / norm);
        }

        var basis = ComplexMatrix.Build.Dense(dim, orthonormal.Count);
        for (int i = 0; i < orthonormal.Count; i++) basis.SetColumn(i, orthonormal[i]);

        // Verify orthonormality after Gram-Schmidt: gram = B^† B should be identity.
        var gram = basis.ConjugateTranspose() * basis;
        double residual = 0;
        for (int i = 0; i < gram.RowCount; i++)
            for (int j = 0; j < gram.ColumnCount; j++)
            {
                double target = i == j ? 1.0 : 0.0;
                double mag = (gram[i, j] - target).Magnitude;
                if (mag > residual) residual = mag;
            }

        return new MultiKBasis(block, quartets, basis, residual);
    }

    /// <summary>Project a full-block matrix onto the multi-k basis: B† · M · B.</summary>
    public ComplexMatrix Project(ComplexMatrix m) =>
        BasisMatrix.ConjugateTranspose() * m * BasisMatrix;

    /// <summary>Project a full-block vector onto the multi-k basis: B† · v.</summary>
    public ComplexVector Project(ComplexVector v) =>
        BasisMatrix.ConjugateTranspose() * v;
}

/// <summary>Structural metadata for the k-th quartet: HD pair (2k−1, 2k+1) coupling and
/// its top singular value σ_0^{(k)}. The pre-Gram-Schmidt vectors are not preserved —
/// after orthonormalisation the basis columns are linear combinations of the original
/// quartet vectors, so the structural identifiers (K, Hd1, Hd2, Sigma0) are the only
/// post-GS-meaningful fields.</summary>
public sealed record KQuartet(int K, int Hd1, int Hd2, double Sigma0);
