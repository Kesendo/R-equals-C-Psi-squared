using System.Collections.Generic;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The site-summed spectator W(ρ) = Σ_l c_l†ρc_l as a matrix between joint-popcount blocks: the
/// shared builder behind Theorem B of <c>docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md</c> (typed:
/// <c>SpectatorIntertwinerClaim</c>). Consumed by BOTH the gate
/// (<c>SpectatorIntertwinerGateTests</c>, SLOW_MSM) and the live witness (<see cref="SectorBraidWitness"/>
/// node 2), so the gate's machine-zero numbers and the inspect-time numbers come from one construction.
///
/// <para><b>Conventions</b> (identical to the gate / <c>SectorBraidModeGeometry.BuildBlock</c>): block basis =
/// <see cref="BlockBasis.PopcountStates"/> ascending, flat index = pIdx·Mq + qIdx. Site indexing is big-endian
/// (site 0 = MSB, site l = bit N−1−l), matching <c>JwBlockBasis</c>/<c>XyJordanWignerModes</c>. The
/// Jordan-Wigner string sign is s_l(a) = (−1)^(# occupied sites m &lt; l in a), i.e. the parity of the bits
/// strictly above bit N−1−l (<see cref="JwSign"/>). On a coherence, W(|a⟩⟨b|) = Σ_{l∉a ∧ l∉b}
/// s_l(a)·s_l(b)·|a∪e_l⟩⟨b∪e_l|: the SAME site l is added to bra and ket, which is exactly why the dephasing
/// sign squares away (Lemma 2) while the single-mode spectator's cross-site signs do not.</para></summary>
public static class SpectatorIntertwiner
{
    /// <summary>Site <paramref name="site"/> occupied in computational-basis state (big-endian: site 0 = MSB
    /// = bit n−1).</summary>
    public static bool Occupied(int n, long state, int site) => ((state >> (n - 1 - site)) & 1) != 0;

    /// <summary>The single-site bit mask of <paramref name="site"/> (big-endian).</summary>
    public static long SiteMask(int n, int site) => 1L << (n - 1 - site);

    /// <summary>JW string sign s_l(a) = (−1)^(# occupied sites strictly before site l), sites ordered 0..n−1
    /// with site 0 = MSB: the sites m &lt; l occupy the bits strictly above bit n−1−l, i.e. state &gt;&gt; (n − l).</summary>
    public static int JwSign(int n, long state, int site)
        => (BitOperations.PopCount((ulong)(state >> (n - site))) & 1) == 0 ? 1 : -1;

    /// <summary>The site-summed spectator W as a matrix from block (pIn, qIn) to block (pIn+1, qIn+1),
    /// in the <see cref="BlockBasis.PopcountStates"/> ordering (flat = pIdx·Mq + qIdx). For each input
    /// coherence |a⟩⟨b| and each site l with l∉a AND l∉b: coefficient s_l(a)·s_l(b) into |a∪e_l⟩⟨b∪e_l|
    /// (the bra-side JW sign is real, so no conjugation enters).</summary>
    public static Matrix<Complex> BuildW(int n, int pIn, int qIn)
    {
        var inP = BlockBasis.PopcountStates(n, pIn);
        var inQ = BlockBasis.PopcountStates(n, qIn);
        var outP = BlockBasis.PopcountStates(n, pIn + 1);
        var outQ = BlockBasis.PopcountStates(n, qIn + 1);
        var outPIdx = new Dictionary<long, int>();
        for (int i = 0; i < outP.Count; i++) outPIdx[outP[i]] = i;
        var outQIdx = new Dictionary<long, int>();
        for (int i = 0; i < outQ.Count; i++) outQIdx[outQ[i]] = i;

        var w = Matrix<Complex>.Build.Dense(outP.Count * outQ.Count, inP.Count * inQ.Count);
        for (int i = 0; i < inP.Count; i++)
            for (int j = 0; j < inQ.Count; j++)
            {
                long a = inP[i], b = inQ[j];
                int col = i * inQ.Count + j;
                for (int l = 0; l < n; l++)
                {
                    if (Occupied(n, a, l) || Occupied(n, b, l)) continue;
                    int row = outPIdx[a | SiteMask(n, l)] * outQ.Count + outQIdx[b | SiteMask(n, l)];
                    w[row, col] += JwSign(n, a, l) * JwSign(n, b, l);
                }
            }
        return w;
    }
}
