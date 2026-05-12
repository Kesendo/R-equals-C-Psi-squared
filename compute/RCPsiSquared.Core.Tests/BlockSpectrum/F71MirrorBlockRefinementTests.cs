using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Tests.TestHelpers;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

public class F71MirrorBlockRefinementTests
{
    // ----------------------------------------------------------------------
    // Claim metadata
    // ----------------------------------------------------------------------

    [Fact]
    public void F71MirrorBlockRefinement_Tier_IsTier1Derived()
    {
        var claim = new F71MirrorBlockRefinement(new JointPopcountSectors());
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void F71MirrorBlockRefinement_DisplayName_IsNonEmpty()
    {
        var claim = new F71MirrorBlockRefinement(new JointPopcountSectors());
        Assert.False(string.IsNullOrWhiteSpace(claim.DisplayName));
    }

    [Fact]
    public void F71MirrorBlockRefinement_Summary_IsNonEmpty()
    {
        var claim = new F71MirrorBlockRefinement(new JointPopcountSectors());
        Assert.False(string.IsNullOrWhiteSpace(claim.Summary));
    }

    [Fact]
    public void F71MirrorBlockRefinement_Constructor_NullParent_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => new F71MirrorBlockRefinement(null!));
    }

    // ----------------------------------------------------------------------
    // RefineWithF71: structural invariants
    // ----------------------------------------------------------------------

    [Fact]
    public void F71MirrorBlockRefinement_RefineWithF71_NullDecomp_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => F71MirrorBlockRefinement.RefineWithF71(null!));
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void F71MirrorBlockRefinement_BasisChange_IsOrthogonal(int N)
    {
        var baseDecomp = JointPopcountSectorBuilder.Build(N);
        var refined = F71MirrorBlockRefinement.RefineWithF71(baseDecomp);

        // Q^T · Q = I (real orthogonal); using ConjugateTranspose for the complex matrix.
        var qtq = refined.BasisChange.ConjugateTranspose() * refined.BasisChange;
        int dim = refined.D * refined.D;
        var ident = Matrix<Complex>.Build.DenseIdentity(dim);
        double err = (qtq - ident).FrobeniusNorm();
        Assert.True(err < 1e-10,
            $"N={N}: Q^T·Q deviation from identity = {err:E3} (expected real-orthogonal Q).");
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void F71MirrorBlockRefinement_SubBlockSizes_SumTo4ToTheN(int N)
    {
        var baseDecomp = JointPopcountSectorBuilder.Build(N);
        var refined = F71MirrorBlockRefinement.RefineWithF71(baseDecomp);
        long total = refined.SectorRanges.Sum(s => (long)s.Size);
        Assert.Equal(1L << (2 * N), total);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void F71MirrorBlockRefinement_SubBlockOffsets_ContiguousNonEmpty(int N)
    {
        // Non-empty sub-blocks must occupy contiguous, non-overlapping column ranges
        // of refined.BasisChange.
        var baseDecomp = JointPopcountSectorBuilder.Build(N);
        var refined = F71MirrorBlockRefinement.RefineWithF71(baseDecomp);

        int expectedOffset = 0;
        foreach (var sub in refined.SectorRanges)
        {
            if (sub.Size == 0) continue;
            Assert.Equal(expectedOffset, sub.Offset);
            expectedOffset += sub.Size;
        }
        Assert.Equal(refined.D * refined.D, expectedOffset);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void F71MirrorBlockRefinement_SectorCount_Is2xJointPopcount(int N)
    {
        // Each joint-popcount sector emits one even + one odd entry (some odd entries may
        // be empty when a sector is entirely F71-fixed); total entries = 2 · (N+1)².
        var baseDecomp = JointPopcountSectorBuilder.Build(N);
        var refined = F71MirrorBlockRefinement.RefineWithF71(baseDecomp);
        Assert.Equal(2 * JointPopcountSectors.SectorCount(N), refined.SectorRanges.Count);
    }

    // ----------------------------------------------------------------------
    // Bit-exact off-block Frobenius after F71 refinement
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void F71MirrorBlockRefinement_OffBlockFrobenius_IsZero_ChainXYZDeph(int N)
    {
        const double J = 1.0;
        const double gamma = 0.5;
        var L = BuildXYZDephasingL(N, J, gamma);

        var baseDecomp = JointPopcountSectorBuilder.Build(N);
        var refined = F71MirrorBlockRefinement.RefineWithF71(baseDecomp);

        // Build the full transformed matrix Lp = Q^T L Q, then verify off-sub-block entries are 0.
        var Lp = refined.BasisChange.ConjugateTranspose() * L * refined.BasisChange;

        double offBlockFroSq = 0.0;
        var subs = refined.SectorRanges;
        for (int si = 0; si < subs.Count; si++)
        {
            for (int sj = 0; sj < subs.Count; sj++)
            {
                if (si == sj) continue;
                if (subs[si].Size == 0 || subs[sj].Size == 0) continue;
                for (int rOff = 0; rOff < subs[si].Size; rOff++)
                {
                    int rowFlat = subs[si].Offset + rOff;
                    for (int cOff = 0; cOff < subs[sj].Size; cOff++)
                    {
                        int colFlat = subs[sj].Offset + cOff;
                        var z = Lp[rowFlat, colFlat];
                        offBlockFroSq += z.Real * z.Real + z.Imaginary * z.Imaginary;
                    }
                }
            }
        }
        double offBlockFro = Math.Sqrt(offBlockFroSq);
        Assert.True(offBlockFro < 1e-10,
            $"N={N}: F71-refined off-sub-block Frobenius = {offBlockFro:E3} (expected ~0; chain L commutes with P_F71⊗P_F71).");
    }

    // ----------------------------------------------------------------------
    // Spectrum bit-exact match: F71-refined per-block eig vs full-L eig
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void F71MirrorBlockRefinement_ComputeSpectrum_MatchesFullL_Multiset(int N)
    {
        var L = BuildXYZDephasingL(N, J: 1.0, gamma: 0.5);

        var eigsFull = L.Evd().EigenValues.ToArray();
        var eigsRefined = F71MirrorBlockRefinement.ComputeSpectrum(L, N);

        Assert.Equal(eigsFull.Length, eigsRefined.Length);
        MultisetAssert.NearestNeighbourEqual(eigsFull, eigsRefined, tolerance: 1e-9, context: $"N={N}");
    }

    [Fact]
    public void F71MirrorBlockRefinement_ComputeSpectrum_NullL_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => F71MirrorBlockRefinement.ComputeSpectrum(null!, 3));
    }

    [Fact]
    public void F71MirrorBlockRefinement_ComputeSpectrum_WrongDimension_Throws()
    {
        var wrong = Matrix<Complex>.Build.DenseIdentity(16);
        Assert.Throws<ArgumentException>(() => F71MirrorBlockRefinement.ComputeSpectrum(wrong, 3));
    }

    // ----------------------------------------------------------------------
    // N=8 projection: max sub-block size approximately MaxSectorSize(8) / 2
    // ----------------------------------------------------------------------

    [Fact]
    public void F71MirrorBlockRefinement_N8_MaxSubBlockSize_ApproxHalvesMaxSectorSize()
    {
        // We don't actually diagonalise at N=8, just count orbit sizes inside the (4, 4)
        // joint-popcount sector (the largest one, size 4900). The F71-even half = (#fixed
        // pairs + #non-fixed pairs); F71-odd half = #non-fixed pairs. At N=8 the central
        // (4, 4) sector has both col and row drawn from C(8, 4) = 70 binary patterns.
        const int N = 8;
        const int d = 1 << N;
        const int mid = N / 2;

        // Enumerate all (col, row) pairs in the (4, 4) sector, counting F71-fixed ones
        // (col palindromic AND row palindromic) and non-fixed ones (paired by the F71
        // involution).
        int fixedCount = 0;
        int totalCount = 0;
        var paired = new HashSet<int>();
        var palindromic = new bool[d];
        for (int x = 0; x < d; x++)
        {
            int m = 0;
            for (int i = 0; i < N; i++)
                if (((x >> i) & 1) != 0) m |= 1 << (N - 1 - i);
            palindromic[x] = (m == x);
        }
        for (int row = 0; row < d; row++)
        {
            if (System.Numerics.BitOperations.PopCount((uint)row) != mid) continue;
            for (int col = 0; col < d; col++)
            {
                if (System.Numerics.BitOperations.PopCount((uint)col) != mid) continue;
                totalCount++;
                if (palindromic[row] && palindromic[col]) fixedCount++;
            }
        }
        Assert.Equal(4900, totalCount);

        int nonFixed = totalCount - fixedCount;
        // nonFixed must be even (paired by F71 involution).
        Assert.Equal(0, nonFixed % 2);
        int evenHalf = fixedCount + nonFixed / 2;
        int oddHalf = nonFixed / 2;

        // Max sub-block at (4, 4) is the F71-even half. Verify it is approximately
        // 4900 / 2 = 2450, with deviation = fixedCount / 2 (the boundary effect).
        Assert.Equal(4900, evenHalf + oddHalf);
        Assert.True(evenHalf <= 4900 / 2 + fixedCount,
            $"F71-even half of (4,4) sector = {evenHalf}; expected ≈ 2450 (= 4900/2); fixedCount = {fixedCount}.");
        Assert.True(evenHalf >= (4900 - fixedCount) / 2,
            $"F71-even half of (4,4) sector = {evenHalf} appears below the lower bound (4900-fixed)/2 = {(4900 - fixedCount) / 2}.");
        // The boundary effect (fixed-point bias) must be small relative to the full sector.
        Assert.True(fixedCount * 4 < 4900,
            $"F71-fixed-point count in (4,4) sector = {fixedCount}; expected ≪ 4900 (boundary effect).");
    }

    // ----------------------------------------------------------------------
    // Per-block path: no full-L materialisation (Phase 4 N=7,8 enabler)
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void F71MirrorBlockRefinement_ComputeSpectrumPerBlock_MatchesFullL(int N)
    {
        var L = BuildXYZDephasingL(N, J: 1.0, gamma: 0.5);
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(0.5, N).ToArray();

        var eigsFull = L.Evd().EigenValues.ToArray();
        var eigsPerBlock = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, gammaPerSite, N);

        Assert.Equal(eigsFull.Length, eigsPerBlock.Length);
        MultisetAssert.NearestNeighbourEqual(eigsFull, eigsPerBlock, tolerance: 1e-9, context: $"N={N}");
    }

    [Fact]
    public void F71MirrorBlockRefinement_ComputeSpectrumPerBlock_NullH_Throws()
    {
        var gammas = new double[3];
        Assert.Throws<ArgumentNullException>(() =>
            F71MirrorBlockRefinement.ComputeSpectrumPerBlock(null!, gammas, 3));
    }

    [Fact]
    public void F71MirrorBlockRefinement_ComputeSpectrumPerBlock_WrongHDim_Throws()
    {
        var wrong = Matrix<Complex>.Build.DenseIdentity(16);
        var gammas = new double[3];
        Assert.Throws<ArgumentException>(() =>
            F71MirrorBlockRefinement.ComputeSpectrumPerBlock(wrong, gammas, 3));
    }

    // ----------------------------------------------------------------------
    // Helpers
    // ----------------------------------------------------------------------

    private static ComplexMatrix BuildXYZDephasingL(int N, double J, double gamma)
    {
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();
        return PauliDephasingDissipator.BuildZ(H, gammaPerSite);
    }
}
