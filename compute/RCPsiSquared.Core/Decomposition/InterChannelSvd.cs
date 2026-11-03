using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Decomposition;

/// <summary>SVD of the inter-channel coupling block of M_H_total.
///
/// <para>For HD₁ = 2k−1 and HD₂ = 2k+1 channels of a coherence block, V_inter =
/// P_{HD₁}^† · M_H_total · P_{HD₂} is the H matrix elements between the two HD subspaces.
/// Its top singular vectors |u_0⟩ ∈ HD₁ and |v_0⟩ ∈ HD₂ (lifted to the full block) are the
/// modes that maximally couple under H — the heuristic "EP-partner" pair from the F86
/// 4-mode picture.</para>
///
/// <para>Empirical findings (c=2 chains N=5..8): σ_0 → 2√2 asymptotically; the probe (Dicke
/// state) is orthogonal to |u_0⟩, |v_0⟩; the 2-level same-sign-imaginary off-diagonal form
/// of the proof is realised here as a ±σ_0 real off-diagonal pattern (different phase
/// convention).</para>
/// </summary>
public sealed class InterChannelSvd
{
    public CoherenceBlock Block { get; }
    public int HammingDistance1 { get; }
    public int HammingDistance2 { get; }
    public IReadOnlyList<double> SingularValues { get; }
    public ComplexVector U0InFullBlock { get; }
    public ComplexVector V0InFullBlock { get; }

    public double Sigma0 => SingularValues[0];

    private InterChannelSvd(CoherenceBlock block, int hd1, int hd2,
        double[] singular, ComplexVector u0, ComplexVector v0)
    {
        Block = block;
        HammingDistance1 = hd1;
        HammingDistance2 = hd2;
        SingularValues = singular;
        U0InFullBlock = u0;
        V0InFullBlock = v0;
    }

    public static InterChannelSvd Build(CoherenceBlock block, int hd1, int hd2)
    {
        var pHd1 = HdSubspaceProjector.Build(block, hd1);
        var pHd2 = HdSubspaceProjector.Build(block, hd2);

        var vInter = pHd1.ConjugateTranspose() * block.Decomposition.MhTotal * pHd2; // (n_hd1 × n_hd2)

        var svd = vInter.Svd();
        var singular = svd.S.Select(z => z.Real).ToArray();
        var (u0Full, v0Full) = TopPair(block, pHd1, pHd2, vInter, svd);
        return new InterChannelSvd(block, hd1, hd2, singular, u0Full, v0Full);
    }

    /// <summary>Relative gap below which a singular value counts as degenerate with σ_0 (the
    /// tolerance of <see cref="F86.Item1Derivation.C2InterChannelProjector.DegeneracyTolerance"/>).</summary>
    public const double DegeneracyTolerance = 1e-8;

    /// <summary>The top inter-channel pair (u_0, v_0). Where σ_0 is simple (odd N at c = 2) this is
    /// the SVD's own top pair. Where it is degenerate (even N at c = 2, twofold; see
    /// <see cref="F86.Item1Derivation.C2InterChannelProjector"/>) a single SVD column is whatever
    /// mixture the library returns, and the four-mode K curve depends on that mixture. The chain
    /// mirror R splits the degenerate space into an R-even and an R-odd direction; the Dicke probe,
    /// the spatial-sum kernel and M_H_total are all R-symmetric, so ρ(t) stays R-even and the R-even
    /// direction is the one the full dynamics uses. u_0 is that direction, v_0 = V†u_0/σ_0 its partner.</summary>
    private static (ComplexVector U0, ComplexVector V0) TopPair(CoherenceBlock block,
        ComplexMatrix pHd1, ComplexMatrix pHd2, ComplexMatrix vInter,
        MathNet.Numerics.LinearAlgebra.Factorization.Svd<Complex> svd)
    {
        double s0 = svd.S[0].Real;
        int deg = 1;
        while (deg < svd.S.Count && Math.Abs(svd.S[deg].Real - s0) / s0 < DegeneracyTolerance) deg++;

        var uRed = svd.U.Column(0);
        if (deg > 1)
        {
            // R-even part of the degenerate space: (u_k + R u_k) for each top column, then its
            // dominant direction. R maps the basis pair (p, q) to (rev p, rev q).
            var basis = block.Basis;
            int n = basis.MTotal;
            var mirror = new int[n];
            foreach (long p in basis.StatesP)
                foreach (long q in basis.StatesQ)
                    mirror[basis.FlatIndex(p, q)] = basis.FlatIndex(Reverse(p, block.N), Reverse(q, block.N));
            var even = ComplexMatrix.Build.Dense(n, deg);
            for (int k = 0; k < deg; k++)
            {
                var u = pHd1 * svd.U.Column(k);
                for (int i = 0; i < n; i++) even[i, k] = u[i] + u[mirror[i]];
            }
            var dominant = even.Svd().U.Column(0);
            uRed = pHd1.ConjugateTranspose() * dominant;
            uRed = uRed / uRed.L2Norm();
        }
        if (deg == 1)
            return (pHd1 * uRed, pHd2 * svd.VT.ConjugateTranspose().Column(0));
        var vRed = vInter.ConjugateTranspose() * uRed;
        return (pHd1 * uRed, pHd2 * (vRed / vRed.L2Norm()));
    }

    private static long Reverse(long state, int n)
    {
        long r = 0;
        for (int i = 0; i < n; i++) if (((state >> i) & 1L) != 0) r |= 1L << (n - 1 - i);
        return r;
    }
}
