using RCPsiSquared.Core.CoherenceBlocks;

namespace RCPsiSquared.Core.F86;

/// <summary>Closed-form C_block(t) for the (popcount-n, popcount-(n+1)) coherence block
/// of an N-qubit chain under uniform-J Heisenberg + local Z-dephasing γ, evaluated on
/// the maximally-coherent initial state |ψ⟩ = (|D_n⟩ + |D_{n+1}⟩)/√2.
///
/// <para>Tier 1 derived (chromaticity-universal):</para>
/// <code>
///     C_block(t)  =  (1/4) · Σ_{k=0..c-1}  (M_{HD=2k+1} / M_block) · exp(−4γ·(2k+1)·t)
/// </code>
///
/// <para>where:
/// <list type="bullet">
///   <item><c>c = min(n, N−1−n) + 1</c> (chromaticity, F74)</item>
///   <item><c>M_block = C(N, n) · C(N, n+1)</c> (block dimension)</item>
///   <item><c>M_{HD=2k+1} = C(N, n−k) · C(N−n+k, k) · C(N−n, k+1)</c> (entries with
///   Hamming distance 2k+1 between popcount-n and popcount-(n+1) states)</item>
/// </list></para>
///
/// <para>At t = 0: <c>Σ M_{HD=2k+1} = M_block</c> by combinatorial identity, giving
/// <c>C_block(0) = M_block/(4·M_block) = 1/4</c> EXACTLY for any (N, c, n) — the
/// algebraic chromaticity-universal anchor on the R=CΨ² Mandelbrot boundary.</para>
///
/// <para>Derivation: max-coherent initial has ρ_ab = 1/(2·√M_block) at every block
/// entry. Pure dephasing decays |ρ_ab(t)|² = |ρ_ab(0)|² · exp(−4γ·HD(a,b)·t). Summing
/// over the block with HD-class counts gives the closed form. The Hamiltonian's
/// contribution is zero on the channel-uniform initial (F73 sum-rule).</para>
/// </summary>
public static class BlockCpsiClosedForm
{
    /// <summary>The Hamming-distance-(2k+1) entry count
    /// <c>C(N, n−k) · C(N−n+k, k) · C(N−n, k+1)</c> for the (popcount-n, popcount-(n+1))
    /// block. Returns 0 when k violates k ≤ n or k+1 ≤ N−n.</summary>
    public static long HammingDistanceCount(int N, int n, int k)
    {
        if (k < 0 || k > n || k + 1 > N - n) return 0;
        return Binomial(N, n - k) * Binomial(N - n + k, k) * Binomial(N - n, k + 1);
    }

    /// <summary>Block dimension M_block = C(N, n) · C(N, n+1).</summary>
    public static long BlockDimension(int N, int n) => Binomial(N, n) * Binomial(N, n + 1);

    /// <summary>Closed-form C_block(t) at the given time. Chromaticity is computed
    /// internally from (N, n) via <see cref="Chromaticity.Compute"/>.</summary>
    public static double At(int N, int n, double gammaZero, double t)
    {
        int c = Chromaticity.Compute(N, n);
        long mBlock = BlockDimension(N, n);
        double sum = 0.0;
        for (int k = 0; k < c; k++)
        {
            long mHd = HammingDistanceCount(N, n, k);
            double weight = (double)mHd / (4.0 * mBlock);
            double rate = 4.0 * gammaZero * (2 * k + 1);
            sum += weight * Math.Exp(-rate * t);
        }
        return sum;
    }

    private static long Binomial(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        if (k == 0 || k == n) return 1;
        if (k > n - k) k = n - k;
        long result = 1;
        for (int i = 0; i < k; i++)
        {
            result *= (n - i);
            result /= (i + 1);
        }
        return result;
    }
}
