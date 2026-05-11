using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

public class LiouvillianBlockSpectrumTests
{
    // ----------------------------------------------------------------------
    // Claim metadata
    // ----------------------------------------------------------------------

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var claim = new LiouvillianBlockSpectrum(new JointPopcountSectors());
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void DisplayName_IsNonEmpty()
    {
        var claim = new LiouvillianBlockSpectrum(new JointPopcountSectors());
        Assert.False(string.IsNullOrWhiteSpace(claim.DisplayName));
    }

    [Fact]
    public void Summary_IsNonEmpty()
    {
        var claim = new LiouvillianBlockSpectrum(new JointPopcountSectors());
        Assert.False(string.IsNullOrWhiteSpace(claim.Summary));
    }

    [Fact]
    public void Constructor_NullParent_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => new LiouvillianBlockSpectrum(null!));
    }

    // ----------------------------------------------------------------------
    // Argument validation
    // ----------------------------------------------------------------------

    [Fact]
    public void ComputeSpectrum_NullL_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => LiouvillianBlockSpectrum.ComputeSpectrum(null!, 3));
    }

    [Fact]
    public void ComputeSpectrum_WrongDimension_Throws()
    {
        // N=3 expects 4^3 = 64; pass a 16×16 (matches N=2) → should throw.
        var wrong = Matrix<Complex>.Build.DenseIdentity(16);
        Assert.Throws<ArgumentException>(() => LiouvillianBlockSpectrum.ComputeSpectrum(wrong, 3));
    }

    // ----------------------------------------------------------------------
    // Spectrum bit-exact match: per-block eig vs direct full-L eig
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void ComputeSpectrum_MatchesFullL_Multiset_AtUnitJAndHalfGamma(int N)
    {
        var L = BuildXYZDephasingL(N, J: 1.0, gamma: 0.5);

        var eigsFull = L.Evd().EigenValues.ToArray();
        var eigsBlock = LiouvillianBlockSpectrum.ComputeSpectrum(L, N);

        Assert.Equal(eigsFull.Length, eigsBlock.Length);
        AssertMultisetEqual(eigsFull, eigsBlock, tol: 1e-9, N: N);
    }

    [Theory]
    [InlineData(3, 0.1)]
    [InlineData(3, 0.5)]
    [InlineData(3, 2.0)]
    public void ComputeSpectrum_MatchesFullL_VaryingGamma_AtN3(int N, double gamma)
    {
        var L = BuildXYZDephasingL(N, J: 1.0, gamma: gamma);

        var eigsFull = L.Evd().EigenValues.ToArray();
        var eigsBlock = LiouvillianBlockSpectrum.ComputeSpectrum(L, N);

        AssertMultisetEqual(eigsFull, eigsBlock, tol: 1e-9, N: N);
    }

    [Theory]
    [InlineData(3, 0.5)]
    [InlineData(3, 1.0)]
    [InlineData(3, 3.0)]
    public void ComputeSpectrum_MatchesFullL_VaryingJ_AtN3(int N, double J)
    {
        var L = BuildXYZDephasingL(N, J: J, gamma: 0.5);

        var eigsFull = L.Evd().EigenValues.ToArray();
        var eigsBlock = LiouvillianBlockSpectrum.ComputeSpectrum(L, N);

        AssertMultisetEqual(eigsFull, eigsBlock, tol: 1e-9, N: N);
    }

    // ----------------------------------------------------------------------
    // Structural checks: sector count and block sizing
    // ----------------------------------------------------------------------

    [Fact]
    public void ComputeSpectrum_AtN4_Produces256Eigenvalues_Across25Blocks()
    {
        const int N = 4;
        var L = BuildXYZDephasingL(N, J: 1.0, gamma: 0.5);

        var eigsBlock = LiouvillianBlockSpectrum.ComputeSpectrum(L, N);

        // 4^4 = 256 total eigenvalues
        Assert.Equal(256, eigsBlock.Length);

        // (N+1)^2 = 25 sectors
        Assert.Equal(25, JointPopcountSectors.SectorCount(N));

        // Total sector size also sums to 256
        var decomp = JointPopcountSectorBuilder.Build(N);
        int totalSize = decomp.SectorRanges.Sum(s => s.Size);
        Assert.Equal(256, totalSize);
        Assert.Equal(25, decomp.SectorRanges.Count);
    }

    [Fact]
    public void ComputeSpectrum_AtN5_LargestBlockIs100()
    {
        const int N = 5;
        var decomp = JointPopcountSectorBuilder.Build(N);
        int maxBlock = decomp.SectorRanges.Max(s => s.Size);

        // C(5, 2)^2 = 10^2 = 100
        Assert.Equal(100, maxBlock);
        Assert.Equal(100L, JointPopcountSectors.MaxSectorSize(N));
    }

    // ----------------------------------------------------------------------
    // Per-block path: no full-L materialisation (Phase 4 N=7,8 enabler)
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void ComputeSpectrumPerBlock_MatchesFullL_AtUnitJAndHalfGamma(int N)
    {
        var L = BuildXYZDephasingL(N, J: 1.0, gamma: 0.5);
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(0.5, N).ToArray();

        var eigsFull = L.Evd().EigenValues.ToArray();
        var eigsPerBlock = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N);

        Assert.Equal(eigsFull.Length, eigsPerBlock.Length);
        AssertMultisetEqual(eigsFull, eigsPerBlock, tol: 1e-9, N: N);
    }

    [Fact]
    public void ComputeSpectrumPerBlock_NullH_Throws()
    {
        var gammas = new double[3];
        Assert.Throws<ArgumentNullException>(() =>
            LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(null!, gammas, 3));
    }

    [Fact]
    public void ComputeSpectrumPerBlock_WrongHDim_Throws()
    {
        // N=3 expects 2^3 = 8; pass a 16×16 → mismatch.
        var wrong = Matrix<Complex>.Build.DenseIdentity(16);
        var gammas = new double[3];
        Assert.Throws<ArgumentException>(() =>
            LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(wrong, gammas, 3));
    }

    [Fact]
    public void ComputeSpectrumPerBlock_WrongGammaLength_Throws()
    {
        const int N = 3;
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var gammas = new double[N + 1];  // wrong length
        Assert.Throws<ArgumentException>(() =>
            LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammas, N));
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

    /// <summary>Assert that two complex-valued multisets are equal up to <paramref name="tol"/>.
    /// Uses greedy nearest-neighbour matching rather than sort-then-compare: Liouvillian
    /// eigenvalues come in highly-degenerate Re-tied clusters (especially under JW for the
    /// XY chain), so a deterministic key sort over (Re, Im) can swap conjugate-pair members
    /// between the two paths even when the multisets are bit-identical.
    /// For each expected eigenvalue, finds the nearest still-unmatched actual eigenvalue
    /// and pairs them, asserting the distance is within tolerance.</summary>
    private static void AssertMultisetEqual(IReadOnlyList<Complex> expected, IReadOnlyList<Complex> actual, double tol, int N)
    {
        Assert.Equal(expected.Count, actual.Count);
        int n = expected.Count;
        var taken = new bool[n];
        for (int i = 0; i < n; i++)
        {
            double bestDist = double.PositiveInfinity;
            int bestJ = -1;
            for (int j = 0; j < n; j++)
            {
                if (taken[j]) continue;
                double dist = (expected[i] - actual[j]).Magnitude;
                if (dist < bestDist)
                {
                    bestDist = dist;
                    bestJ = j;
                }
            }
            Assert.True(bestJ >= 0 && bestDist < tol,
                $"N={N}: eigenvalue {i} mismatch — expected={expected[i]}, nearest actual={(bestJ >= 0 ? actual[bestJ].ToString() : "<none>")}, |Δ|={bestDist:E3} (tol={tol:E1}).");
            taken[bestJ] = true;
        }
    }
}
