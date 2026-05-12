using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

/// <summary>Phase D informational tests pinning the boundary of <see cref="LiouvillianBlockSpectrum"/>
/// under T1 amplitude damping (Lindblad operator <c>c_l = √γ1 · σ⁻_l</c> per site).
///
/// <para><b>Empirical finding (this file's purpose).</b> T1 is <em>not</em> block-diagonal under
/// the <see cref="JointPopcountSectorBuilder"/> permutation that <see cref="LiouvillianBlockSpectrum"/>
/// uses for chain XY + Z-dephasing. But it is also <em>not</em> spectrum-breaking: T1 is
/// block <em>upper-triangular</em> in the joint-popcount sector ordering, with 0% strictly
/// below the diagonal blocks and a non-trivial fraction (≈11.76% at γ1 = J) strictly above.
/// Eigenvalues of a block-triangular matrix equal the multiset union of its diagonal-block
/// eigenvalues, so <see cref="LiouvillianBlockSpectrum.ComputeSpectrum"/> still returns the
/// correct full-L spectrum bit-exactly for T1.</para>
///
/// <para><b>Math.</b> σ⁻ = (X − iY)/2 = [[0, 1], [0, 0]] flips |1⟩ → |0⟩, lowering Hilbert-side
/// popcount by 1. In the Lindblad super-operator <c>c ρ c† − ½ {c†c, ρ}</c>, the only off-block
/// term is <c>c ⊗ c*</c> (the "jump"), and for σ⁻ this lowers BOTH ket popcount p_c AND bra
/// popcount p_r by 1. The recycle <c>c†c ⊗ I</c> and <c>I ⊗ (c†c)^T</c> are diagonal in the
/// computational basis and stay on-block. So in the
/// <see cref="JointPopcountSectorBuilder"/> ordering (sectors lex-sorted by (p_c, p_r)),
/// every off-block entry connects sector index s_i to s_j with s_i &lt; s_j (the lowered
/// (p_c−1, p_r−1) sector sits earlier in the lex order).</para>
///
/// <para><b>Consequence for callers.</b> <see cref="LiouvillianBlockSpectrum.ComputeSpectrum"/>
/// is <em>safe to use</em> for T1 spectra. But the joint-popcount basis no longer block-diagonalises
/// L: <em>eigenvectors</em>, partial traces of resolvents, channel inversions, and the
/// per-block direct-sum decomposition itself break under T1. Use full-L eig
/// (<c>RCPsiSquared.Compute.Liouvillian</c>) when right- or left-eigenvectors of L are needed
/// under amplitude damping.</para>
///
/// <para>Phase A (commit ded3ca8) anchors the bit-exact equivalence of <see cref="LiouvillianBlockSpectrum"/>
/// vs full-L diagonalisation under XY+Z-dephasing. Phase D pins the surprise: T1 leaks
/// block-triangularly, not block-densely, so the spectrum survives even though the
/// direct-sum structure does not.</para></summary>
public class T1BreakingInformationalTests
{
    private const int N = 4;
    private const double J = 1.0;
    private const double Gamma1 = 1.0;  // chosen so γ1 ~ J → off-block fraction ≈ 11.76%

    [Fact]
    public void T1Dissipator_Breaks_JointPopcount_BlockDiagonality_ButOnlyUpperTriangularly()
    {
        var L = BuildXYPlusPureT1L(N, J, Gamma1);
        var decomp = JointPopcountSectorBuilder.Build(N);
        int liouvilleDim = 1 << (2 * N);
        Assert.Equal(liouvilleDim, L.RowCount);

        // Map each position in the permuted basis to its sector index (0..(N+1)²−1).
        var sectorIndexByPermPos = new int[liouvilleDim];
        for (int s = 0; s < decomp.SectorRanges.Count; s++)
        {
            var sec = decomp.SectorRanges[s];
            for (int k = 0; k < sec.Size; k++)
                sectorIndexByPermPos[sec.Offset + k] = s;
        }

        // Walk the permuted L and split Frobenius² into on-block, strictly-upper-block,
        // and strictly-lower-block contributions.
        double onBlock2 = 0, offUpper2 = 0, offLower2 = 0;
        for (int i = 0; i < liouvilleDim; i++)
        {
            int permI = decomp.Permutation[i];
            int si = sectorIndexByPermPos[i];
            for (int j = 0; j < liouvilleDim; j++)
            {
                int permJ = decomp.Permutation[j];
                int sj = sectorIndexByPermPos[j];
                double m2 = L[permI, permJ].Magnitude * L[permI, permJ].Magnitude;
                if (si == sj) onBlock2 += m2;
                else if (si < sj) offUpper2 += m2;
                else offLower2 += m2;
            }
        }

        double total = onBlock2 + offUpper2 + offLower2;
        double offUpperFraction = offUpper2 / total;
        double offLowerFraction = offLower2 / total;

        // (1) T1 leaks a non-trivial fraction of L into off-block — joint-popcount is NOT
        // block-diagonal under amplitude damping. With γ1 = J = 1 at N=4 this is ≈11.76%.
        Assert.True(offUpperFraction > 0.05,
            $"Expected T1 to leak >5% of L into strictly-upper off-block; got {offUpperFraction:P3}.");

        // (2) But all of the leak sits strictly above the diagonal blocks: T1 is block
        // upper-triangular, never block-dense. The strict-lower fraction must be machine zero.
        Assert.True(offLowerFraction < 1e-15,
            $"Expected zero strict-lower off-block content; got {offLowerFraction:E3}. " +
            $"If non-zero, T1 is no longer purely block-triangular and BlockSpectrum.ComputeSpectrum " +
            $"would no longer return the correct spectrum.");
    }

    [Fact]
    public void BlockSpectrum_OnT1Liouvillian_StillMatches_FullEig_BecauseBlockTriangular()
    {
        var L = BuildXYPlusPureT1L(N, J, Gamma1);

        // Pipeline 1: BlockSpectrum (assumes block diagonal; in fact L is block-triangular,
        // so the per-block eigenvalues still concatenate to the full spectrum as a multiset).
        var blockSpectrum = LiouvillianBlockSpectrum.ComputeSpectrum(L, N);

        // Pipeline 2: full-L direct eig (ground truth).
        var fullSpectrum = L.Evd().EigenValues.ToArray();

        Assert.Equal(blockSpectrum.Length, fullSpectrum.Length);

        // Greedy multiset comparison: max nearest-neighbour distance after pairing.
        var (matched, maxDelta) = TryMatchMultiset(blockSpectrum, fullSpectrum, tolerance: 1e-9);

        // The Phase D surprise: BlockSpectrum still works for SPECTRA under T1 because
        // L is block triangular in joint-popcount, not block dense. Eigenvalues of a
        // block-triangular matrix equal the multiset union of its diagonal-block
        // eigenvalues. Eigenvectors and resolvents do break, but the spectrum survives.
        Assert.True(matched,
            $"Expected BlockSpectrum to still return the correct multiset for T1-Liouvillians " +
            $"because they are block upper-triangular in joint-popcount. Got max |Δλ| = {maxDelta:E3}. " +
            $"If this fails, T1 has filled below-diagonal blocks and BlockSpectrum.ComputeSpectrum " +
            $"can no longer be applied to T1-Liouvillians without re-deriving.");
    }

    // ----------------------------------------------------------------------
    // Helpers
    // ----------------------------------------------------------------------

    /// <summary>Build XY chain + pure T1 amplitude damping per site (no Z-dephasing).
    /// Uses <see cref="T1Dissipator.Build"/> with <c>γ_Z = 0</c> on every site.</summary>
    private static ComplexMatrix BuildXYPlusPureT1L(int N, double J, double gamma1)
    {
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaZ = Enumerable.Repeat(0.0, N).ToArray();
        var gammaT1 = Enumerable.Repeat(gamma1, N).ToArray();
        return T1Dissipator.Build(H, gammaZ, gammaT1);
    }

    /// <summary>Greedy nearest-neighbour multiset comparison. Returns
    /// <c>(matchedWithinTolerance, maxNearestDistance)</c>.</summary>
    private static (bool Matched, double MaxDelta) TryMatchMultiset(
        IReadOnlyList<Complex> a, IReadOnlyList<Complex> b, double tolerance)
    {
        if (a.Count != b.Count) return (false, double.PositiveInfinity);
        int n = a.Count;
        var taken = new bool[n];
        double maxDelta = 0;
        for (int i = 0; i < n; i++)
        {
            double bestDist = double.PositiveInfinity;
            int bestJ = -1;
            for (int j = 0; j < n; j++)
            {
                if (taken[j]) continue;
                double dist = (a[i] - b[j]).Magnitude;
                if (dist < bestDist)
                {
                    bestDist = dist;
                    bestJ = j;
                }
            }
            if (bestJ < 0) return (false, double.PositiveInfinity);
            taken[bestJ] = true;
            if (bestDist > maxDelta) maxDelta = bestDist;
        }
        return (maxDelta < tolerance, maxDelta);
    }
}
