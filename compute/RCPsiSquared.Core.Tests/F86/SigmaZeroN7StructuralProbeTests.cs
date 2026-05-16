using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum.JordanWigner;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.F86.JordanWigner;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>F86 Item 5 structural probe: σ_0(c=2, N=7) = 2√2 bit-exactly (Tier-1-candidate
/// sweet-spot identity per <see cref="RCPsiSquared.Core.F86.SigmaZeroChromaticityScaling"/>).
/// Why? At N=7 the OBC dispersion ε_k = 2J cos(πk/8) has ε_2 = √2 exactly, so 2√2 = 2·ε_2.
/// This probe checks whether the V_inter top singular vectors |u_0⟩, |v_0⟩ have JW-mode
/// content concentrated on mode k=2 (and what about its partner k=N+1−k=6=−ε_2), which would
/// explain the σ_0 = 2·ε_2 identity as a mode-2 pair-creation dominance.</summary>
public class SigmaZeroN7StructuralProbeTests
{
    private readonly ITestOutputHelper _out;
    public SigmaZeroN7StructuralProbeTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Probe_N7_C2_TopSvdVectors_JwModeContent()
    {
        const int N = 7;
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var svd = InterChannelSvd.Build(block, hd1: 1, hd2: 3);
        double exactSweetSpot = 2.0 * Math.Sqrt(2.0);
        _out.WriteLine($"σ_0(c=2, N=7) = {svd.Sigma0:F12}");
        _out.WriteLine($"2√2 exactly     = {exactSweetSpot:F12}");
        _out.WriteLine($"Δ = {Math.Abs(svd.Sigma0 - exactSweetSpot):G3}");

        var modes = XyJordanWignerModes.Build(N, J: 1.0);
        _out.WriteLine($"\nOBC dispersion ε_k = 2 cos(πk/{N + 1}):");
        for (int k = 1; k <= N; k++)
            _out.WriteLine($"  ε_{k} = {modes.Dispersion[k - 1]:F8}");
        _out.WriteLine($"  ε_2 = √2 = {Math.Sqrt(2):F8} → 2·ε_2 = 2√2 = {2 * Math.Sqrt(2):F8}");

        // Project top SVD vectors |u_0⟩ ∈ HD₁ (popcount-1 row, popcount-2 col, HD=1) and
        // |v_0⟩ ∈ HD₂ (HD=3) onto the JW Slater-pair basis for the c=2 sector
        // (p_c=2, p_r=1). The JW basis indexes by mode triples (k, k₁, k₂) with k single,
        // (k₁<k₂) pair: column α of U has slaterRow[ri, k]·slaterCol[ci, (k₁, k₂)].
        var jwBasis = JwSlaterPairBasis.Build(N, pCol: 2, pRow: 1);
        _out.WriteLine($"\nJW Slater-pair basis dim: {jwBasis.U.RowCount} = |HD₁| + |HD₂| (= 7·21 = {N * N * (N - 1) / 2})");

        // Convert the top SVD vectors (which live in the c=2 block of length 7·21 = 147)
        // from computational basis to JW Slater-pair basis: u_0_jw = U^† · u_0_comp.
        // svd.U0InFullBlock is already in the block's basis (dim 147), same indexing as
        // jwBasis.U rows (per JwSlaterPairBasis.FlatIndices).
        var u0Jw = jwBasis.Uinv * svd.U0InFullBlock;
        var v0Jw = jwBasis.Uinv * svd.V0InFullBlock;

        _out.WriteLine("\nTop-amplitude entries of u_0 (HD=1) in JW Slater-pair basis (|u_0_jw[α]|² ≥ 0.01):");
        for (int alpha = 0; alpha < u0Jw.Count; alpha++)
        {
            double mag2 = u0Jw[alpha].Real * u0Jw[alpha].Real + u0Jw[alpha].Imaginary * u0Jw[alpha].Imaginary;
            if (mag2 < 0.01) continue;
            int Lidx = alpha / jwBasis.ColumnModeSubsets.Count;
            int Kidx = alpha % jwBasis.ColumnModeSubsets.Count;
            var L = jwBasis.RowModeSubsets[Lidx];
            var K = jwBasis.ColumnModeSubsets[Kidx];
            _out.WriteLine($"  α={alpha,3}  L={{{string.Join(",", L)}}}, K={{{string.Join(",", K)}}}, " +
                           $"|amplitude|² = {mag2:F4}");
        }

        _out.WriteLine("\nTop-amplitude entries of v_0 (HD=3) in JW Slater-pair basis (|v_0_jw[α]|² ≥ 0.01):");
        for (int alpha = 0; alpha < v0Jw.Count; alpha++)
        {
            double mag2 = v0Jw[alpha].Real * v0Jw[alpha].Real + v0Jw[alpha].Imaginary * v0Jw[alpha].Imaginary;
            if (mag2 < 0.01) continue;
            int Lidx = alpha / jwBasis.ColumnModeSubsets.Count;
            int Kidx = alpha % jwBasis.ColumnModeSubsets.Count;
            var L = jwBasis.RowModeSubsets[Lidx];
            var K = jwBasis.ColumnModeSubsets[Kidx];
            _out.WriteLine($"  α={alpha,3}  L={{{string.Join(",", L)}}}, K={{{string.Join(",", K)}}}, " +
                           $"|amplitude|² = {mag2:F4}");
        }
    }
}
