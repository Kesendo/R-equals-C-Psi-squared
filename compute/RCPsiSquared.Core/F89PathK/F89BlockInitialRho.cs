using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>F89 path-k block-basis helpers and initial-state builder.
///
/// <para>Big-endian basis-state index convention matching Python:
/// <c>|b_0 b_1 ... b_{n-1}⟩ → Σ b_i · 2^(n-1-i)</c>.</para>
///
/// <para><see cref="ComputeRhoBlockZero"/> is the closed-form partial trace of ρ_cc over the
/// (N − n_block) bare sites, returning the block-reduced d × d Hermitian density matrix
/// (d = 2^n_block) at t = 0. Used as the initial condition for path-k block dynamics under
/// <see cref="F89BlockLiouvillian.BuildBlockL"/>.</para>
/// </summary>
public static class F89BlockInitialRho
{
    /// <summary>State-index bit-position weights: site i in an n_block-qubit word contributes
    /// <c>2^(n_block-1-i)</c> to the integer basis index.</summary>
    public static int[] BlockBitPos(int nBlock)
    {
        if (nBlock < 1)
            throw new ArgumentOutOfRangeException(nameof(nBlock), nBlock, "nBlock must be ≥ 1.");
        var weights = new int[nBlock];
        for (int i = 0; i < nBlock; i++) weights[i] = 1 << (nBlock - 1 - i);
        return weights;
    }

    /// <summary>Compute the basis-state integer index from a bit list and (optionally) a
    /// pre-computed weight array from <see cref="BlockBitPos"/>.</summary>
    public static int StateIdx(int[] bits, int[]? bitPos = null)
    {
        bitPos ??= BlockBitPos(bits.Length);
        if (bits.Length != bitPos.Length)
            throw new ArgumentException("bits length must match bitPos length.");
        int idx = 0;
        for (int i = 0; i < bits.Length; i++) idx += bitPos[i] * bits[i];
        return idx;
    }

    /// <summary>ρ_block(0) = Tr_{N-n_block bare}(ρ_cc) for any path-k block at N qubits.
    ///
    /// <para>ρ_cc = (|S_1⟩⟨S_2| + h.c.) / 2 where S_1 is the symmetric single-excitation Dicke
    /// state and S_2 is the symmetric double-excitation Dicke state. Partial trace yields:</para>
    /// <list type="bullet">
    ///   <item>Term 1 (popcount(c) = 0): Σ_{i ∈ block} Σ_{j&lt;k both in block}
    ///         |SE_i^B⟩⟨DE_{jk}^B|</item>
    ///   <item>Term 2 (popcount(c) = 1): N_E · Σ_{j ∈ block} |0^B⟩⟨SE_j^B|, where
    ///         N_E = N − n_block</item>
    /// </list>
    /// <para>Both terms scaled by <c>pre = 1/√(N²(N-1)/2)</c>; final ρ symmetrised via
    /// <c>(ρ + ρ†) / 2</c>.</para>
    /// </summary>
    public static ComplexMatrix ComputeRhoBlockZero(int nBlock, int N)
    {
        if (nBlock < 1)
            throw new ArgumentOutOfRangeException(nameof(nBlock), nBlock, "nBlock must be ≥ 1.");
        if (N < nBlock)
            throw new ArgumentOutOfRangeException(nameof(N), N, $"N must be ≥ nBlock={nBlock}.");

        int d = 1 << nBlock;
        var bitPos = BlockBitPos(nBlock);
        int nE = N - nBlock;
        double pre = 1.0 / Math.Sqrt((double)N * N * (N - 1) / 2.0);

        var rho = Matrix<Complex>.Build.Dense(d, d);

        for (int i = 0; i < nBlock; i++)
        {
            var bits = new int[nBlock];
            bits[i] = 1;
            int idxSe = StateIdx(bits, bitPos);
            for (int j = 0; j < nBlock; j++)
            {
                for (int k = j + 1; k < nBlock; k++)
                {
                    var bitsDe = new int[nBlock];
                    bitsDe[j] = 1;
                    bitsDe[k] = 1;
                    rho[idxSe, StateIdx(bitsDe, bitPos)] += pre;
                }
            }
        }

        for (int j = 0; j < nBlock; j++)
        {
            var bits = new int[nBlock];
            bits[j] = 1;
            rho[0, StateIdx(bits, bitPos)] += pre * nE;
        }

        var rhoH = (rho + rho.ConjugateTranspose()) / 2.0;
        return rhoH;
    }
}
