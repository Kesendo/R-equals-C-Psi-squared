using System.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace RCPsiSquared.Core.F86;

/// <summary>State-level block-purity content C_block on the (popcount-n, popcount-(n+1))
/// coherence block of any density matrix on a 2^N Hilbert space.
///
/// <para>Definition: <c>C_block(ρ, n) = Σ |ρ_{ab}|²</c> with popcount(a) = n,
/// popcount(b) = n+1.</para>
///
/// <para><b>Theorem 2 of PROOF_BLOCK_CPSI_QUARTER (Tier 1 derived):</b>
/// <c>C_block ≤ 1/4</c> for ANY density matrix on the full 2^N Hilbert space, with
/// equality iff ρ is the canonical Dicke superposition (|D_n⟩+|D_{n+1}⟩)/√2 (up to global
/// phase). The 1/4 maxval is the universal Mandelbrot-cardioid ceiling instanced at the
/// c-block level via the bilinear apex argmax/maxval pair (1/2, 1/4) of p·(1−p); see
/// <see cref="Symmetry.QuarterAsBilinearMaxvalClaim"/> and
/// <see cref="Symmetry.ArgmaxMaxvalPairClaim"/>.</para>
///
/// <para>Companion to <see cref="BlockCpsiClosedForm"/>: that helper computes the
/// trajectory C_block(t) on the maximally-coherent Dicke initial state under pure
/// Z-dephasing; this helper takes an arbitrary ρ (e.g. hardware-reconstructed reduced
/// state) and computes its current C_block content.</para>
/// </summary>
public static class BlockCoherenceContent
{
    /// <summary>Block-purity content on the (popcount-n, popcount-(n+1)) block of an
    /// N-qubit density matrix ρ. The matrix dimension must be a power of 2; n must be
    /// in [0, N−1] so that popcount-(n+1) is a valid sector.</summary>
    public static double Compute(Matrix<Complex> rho, int n)
    {
        if (rho is null) throw new ArgumentNullException(nameof(rho));
        int dim = rho.RowCount;
        if (rho.ColumnCount != dim)
            throw new ArgumentException($"rho must be square; got {rho.RowCount}×{rho.ColumnCount}", nameof(rho));
        int N = (int)Math.Round(Math.Log2(dim));
        if ((1 << N) != dim)
            throw new ArgumentException($"rho dimension {dim} is not a power of 2", nameof(rho));
        if (n < 0 || n + 1 > N)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"n must be in [0, {N - 1}] for an N={N} qubit state; got {n}");

        double sum = 0.0;
        for (int a = 0; a < dim; a++)
        {
            if (System.Numerics.BitOperations.PopCount((uint)a) != n) continue;
            for (int b = 0; b < dim; b++)
            {
                if (System.Numerics.BitOperations.PopCount((uint)b) != n + 1) continue;
                double mag = rho[a, b].Magnitude;
                sum += mag * mag;
            }
        }
        return sum;
    }

    /// <summary>The universal Theorem 2 ceiling. C_block ≤ <see cref="Quarter"/> for any
    /// density matrix on 2^N (proof: PROOF_BLOCK_CPSI_QUARTER Theorem 2; argmax/maxval
    /// pair of the bilinear apex p·(1−p), Tier 1 derived).</summary>
    public const double Quarter = 0.25;
}
