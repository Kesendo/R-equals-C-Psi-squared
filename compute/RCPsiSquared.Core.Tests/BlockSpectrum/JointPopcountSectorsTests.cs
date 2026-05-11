using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

public class JointPopcountSectorsTests
{
    // ----------------------------------------------------------------------
    // Claim metadata
    // ----------------------------------------------------------------------

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, new JointPopcountSectors().Tier);
    }

    [Fact]
    public void DisplayName_IsNonEmpty()
    {
        Assert.False(string.IsNullOrWhiteSpace(new JointPopcountSectors().DisplayName));
    }

    [Fact]
    public void Summary_IsNonEmpty()
    {
        Assert.False(string.IsNullOrWhiteSpace(new JointPopcountSectors().Summary));
    }

    // ----------------------------------------------------------------------
    // Combinatorial closed forms
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(1, 4)]    // (1+1)^2 = 4
    [InlineData(2, 9)]    // 3^2 = 9
    [InlineData(3, 16)]
    [InlineData(4, 25)]
    [InlineData(5, 36)]
    [InlineData(6, 49)]
    public void SectorCount_IsNPlus1Squared(int N, int expected)
    {
        Assert.Equal(expected, JointPopcountSectors.SectorCount(N));
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void SectorSizes_SumTo4ToTheN(int N)
    {
        long total = 0;
        for (int pc = 0; pc <= N; pc++)
            for (int pr = 0; pr <= N; pr++)
                total += JointPopcountSectors.SectorSize(N, pc, pr);
        long expected = 1L << (2 * N); // 4^N
        Assert.Equal(expected, total);
    }

    [Theory]
    // C(N, N/2)^2 for N=1..8
    [InlineData(1, 1)]   // C(1,0)^2 = 1; ⌊1/2⌋=0 so C(1,0)=1
    [InlineData(2, 4)]   // C(2,1)^2 = 4
    [InlineData(3, 9)]   // C(3,1)^2 = 9; ⌊3/2⌋=1, C(3,1)=3
    [InlineData(4, 36)]  // C(4,2)^2 = 36
    [InlineData(5, 100)] // C(5,2)^2 = 100
    [InlineData(6, 400)] // C(6,3)^2 = 400
    [InlineData(7, 1225)] // C(7,3)^2 = 35^2 = 1225
    [InlineData(8, 4900)] // C(8,4)^2 = 70^2 = 4900
    public void MaxSectorSize_MatchesBinomialMidSquared(int N, long expected)
    {
        Assert.Equal(expected, JointPopcountSectors.MaxSectorSize(N));
    }

    [Fact]
    public void N8_Projection_Matches4900And81()
    {
        Assert.Equal(4900, JointPopcountSectors.MaxSectorSize(8));
        Assert.Equal(81, JointPopcountSectors.SectorCount(8));
    }

    [Fact]
    public void SectorSize_OutOfRange_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => JointPopcountSectors.SectorSize(3, -1, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => JointPopcountSectors.SectorSize(3, 0, 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => JointPopcountSectors.SectorSize(-1, 0, 0));
    }

    // ----------------------------------------------------------------------
    // Builder structural invariants
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void Build_PermutationIsAValidPermutation(int N)
    {
        var decomp = JointPopcountSectorBuilder.Build(N);
        int liouvilleDim = decomp.D * decomp.D;
        Assert.Equal(liouvilleDim, decomp.Permutation.Length);
        var seen = new bool[liouvilleDim];
        foreach (var idx in decomp.Permutation)
        {
            Assert.InRange(idx, 0, liouvilleDim - 1);
            Assert.False(seen[idx], $"index {idx} appears twice in permutation");
            seen[idx] = true;
        }
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void Build_SectorCountMatchesNPlus1Squared(int N)
    {
        var decomp = JointPopcountSectorBuilder.Build(N);
        Assert.Equal(JointPopcountSectors.SectorCount(N), decomp.SectorRanges.Count);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void Build_SectorSizesMatchClosedForm(int N)
    {
        var decomp = JointPopcountSectorBuilder.Build(N);
        long totalCovered = 0;
        foreach (var sec in decomp.SectorRanges)
        {
            long expected = JointPopcountSectors.SectorSize(N, sec.PCol, sec.PRow);
            Assert.Equal(expected, sec.Size);
            totalCovered += sec.Size;
        }
        Assert.Equal(decomp.D * decomp.D, (int)totalCovered);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void Build_SectorOffsetsAreContiguousAndOrdered(int N)
    {
        var decomp = JointPopcountSectorBuilder.Build(N);
        int expectedOffset = 0;
        foreach (var sec in decomp.SectorRanges)
        {
            Assert.Equal(expectedOffset, sec.Offset);
            expectedOffset += sec.Size;
        }
        Assert.Equal(decomp.D * decomp.D, expectedOffset);
    }

    [Fact]
    public void Build_N_OutOfRange_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => JointPopcountSectorBuilder.Build(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => JointPopcountSectorBuilder.Build(13));
        Assert.Throws<ArgumentOutOfRangeException>(() => JointPopcountSectorBuilder.Build(-1));
    }

    // ----------------------------------------------------------------------
    // Block-diagonal verification (Tier1Derived bit-exact at N=2, 3, 4)
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void XYZDephasing_Liouvillian_IsBlockDiagonalInJointPopcount(int N)
    {
        // Build L for an XY chain (H = (J/2)·Σ_b (X_bX_{b+1}+Y_bY_{b+1})) with uniform
        // per-site Z-dephasing. Convention matches RCPsiSquared.Core.Lindblad
        // (row-major vec(ρ), flat = row*d + col).
        const double J = 1.0;
        const double gamma = 0.5;
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();
        var L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);

        // Build the joint-popcount sector permutation and verify L is exactly
        // block-diagonal once rows/cols are reordered by joint-popcount labels.
        var decomp = JointPopcountSectorBuilder.Build(N);

        // For each (sector_i, sector_j) with i ≠ j (off-diagonal block), the
        // sub-matrix L[perm[i_block], perm[j_block]] must be exactly zero.
        // We accumulate the off-block Frobenius norm.
        double offBlockFroSq = 0.0;
        var sectors = decomp.SectorRanges;
        for (int si = 0; si < sectors.Count; si++)
        {
            for (int sj = 0; sj < sectors.Count; sj++)
            {
                if (si == sj) continue;
                var rs = sectors[si];
                var cs = sectors[sj];
                for (int rOff = 0; rOff < rs.Size; rOff++)
                {
                    int rowFlat = decomp.Permutation[rs.Offset + rOff];
                    for (int cOff = 0; cOff < cs.Size; cOff++)
                    {
                        int colFlat = decomp.Permutation[cs.Offset + cOff];
                        var z = L[rowFlat, colFlat];
                        offBlockFroSq += z.Real * z.Real + z.Imaginary * z.Imaginary;
                    }
                }
            }
        }
        double offBlockFro = Math.Sqrt(offBlockFroSq);
        Assert.True(offBlockFro < 1e-10,
            $"N={N}: off-block Frobenius = {offBlockFro:E3} (expected ~0 for joint-popcount block-diagonality).");
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void XYZDephasing_Liouvillian_TotalFrobeniusSplitsAcrossDiagonalBlocks(int N)
    {
        // Stronger check: Σ ||L_block||_F^2 == ||L||_F^2 to floating-point precision,
        // which together with the off-block test confirms exact block-diagonality.
        const double J = 1.0;
        const double gamma = 0.5;
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();
        var L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);
        var decomp = JointPopcountSectorBuilder.Build(N);

        double fullFroSq = 0.0;
        int liouvilleDim = decomp.D * decomp.D;
        for (int i = 0; i < liouvilleDim; i++)
            for (int j = 0; j < liouvilleDim; j++)
                fullFroSq += L[i, j].Real * L[i, j].Real + L[i, j].Imaginary * L[i, j].Imaginary;

        double diagBlockFroSq = 0.0;
        foreach (var sec in decomp.SectorRanges)
        {
            for (int rOff = 0; rOff < sec.Size; rOff++)
            {
                int rowFlat = decomp.Permutation[sec.Offset + rOff];
                for (int cOff = 0; cOff < sec.Size; cOff++)
                {
                    int colFlat = decomp.Permutation[sec.Offset + cOff];
                    var z = L[rowFlat, colFlat];
                    diagBlockFroSq += z.Real * z.Real + z.Imaginary * z.Imaginary;
                }
            }
        }
        Assert.Equal(fullFroSq, diagBlockFroSq, precision: 10);
    }
}
