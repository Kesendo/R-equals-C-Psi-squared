using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using Xunit;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Guards <see cref="SectorBlock"/>: the RAW joint-popcount (p, q̃) block of the
/// XY (Δ=0) + Z-dephasing Liouvillian L(q)/γ = −i·q·[Ĥ,·] + D, Ĥ = Σ_bonds (XX + YY).
///
/// <para>The RAW block carries NO spatial-symmetry reduction (it is a genuine principal
/// submatrix of the full L in the popcount basis), so it is validated ONLY against the full
/// L decomposition (footing-independent): since L is block-diagonal in the joint-popcount
/// basis, the multiset union of every (p, q̃) block's spectrum must equal the full-L
/// spectrum exactly.</para></summary>
public class SectorBlockTests
{
    // dim of the (p, q̃) block = C(N,p)·C(N,q̃).
    [Theory]
    [InlineData(4, 1, 2, 24)]   // C(4,1)·C(4,2) = 4·6
    [InlineData(4, 0, 1, 4)]    // C(4,0)·C(4,1) = 1·4
    [InlineData(4, 2, 2, 36)]   // C(4,2)·C(4,2) = 6·6
    public void Build_HasExpectedDimension(int N, int p, int qTilde, int expectedDim)
    {
        var block = SectorBlock.Build(N, p, qTilde, new Complex(2.0, 0.0));
        Assert.Equal(expectedDim, block.GetLength(0));
        Assert.Equal(expectedDim, block.GetLength(1));
    }

    // L is block-diagonal in the (p, q̃) basis, so each sector block's eigenvalues are a
    // genuine subset of the full L(2) spectrum, and the union of all (N+1)² blocks is the
    // whole 4^N spectrum. Matched by greedy nearest-neighbour with removal (NOT sort-then-pair:
    // near-degenerate Liouvillian spectra break sorting).
    [Fact]
    public void AllSectorBlocks_UnionToFullLiouvillianSpectrum_At_N4()
    {
        const int N = 4;
        var q = new Complex(2.0, 0.0);

        var full = EigCommon.Eigenvalues(SectorBlock.BuildFullLiouvillian(N, q));
        Assert.Equal(256, full.Length);   // 4^4

        var union = new List<Complex>();
        for (int p = 0; p <= N; p++)
            for (int qTilde = 0; qTilde <= N; qTilde++)
                union.AddRange(EigCommon.Eigenvalues(SectorBlock.Build(N, p, qTilde, q)));

        Assert.Equal(256, union.Count);
        EigCommon.AssertSpectraMatch(full, union.ToArray(), 1e-9);
    }
}

/// <summary>Small eigenvalue helpers for the block-spectrum tests: MathNet dense EVD to a
/// flat <see cref="Complex"/>[] (works for the non-Hermitian Liouvillian, mirroring the
/// trusted <c>block.Evd().EigenValues</c> path in <c>LiouvillianSectorSweep</c>), plus a
/// greedy nearest-neighbour-with-removal multiset comparison robust to near-degeneracy.</summary>
internal static class EigCommon
{
    public static Complex[] Eigenvalues(ComplexMatrix m) => m.Evd().EigenValues.ToArray();

    public static Complex[] Eigenvalues(Complex[,] a) =>
        Eigenvalues(Matrix<Complex>.Build.DenseOfArray(a));

    /// <summary>Assert two spectra are equal as multisets via greedy nearest-neighbour with
    /// removal: for each expected λ, take the closest not-yet-matched actual λ and require it
    /// within <paramref name="tol"/>. NOT Hausdorff (which would pass a spectrum against a
    /// sub-multiset of itself) and NOT sort-then-pair (near-degenerate spectra flip pairings).</summary>
    public static void AssertSpectraMatch(Complex[] expected, Complex[] actual, double tol)
    {
        Assert.Equal(expected.Length, actual.Length);
        var remaining = new List<Complex>(actual);
        foreach (var e in expected)
        {
            int best = -1;
            double bestDist = double.MaxValue;
            for (int j = 0; j < remaining.Count; j++)
            {
                double dist = (remaining[j] - e).Magnitude;
                if (dist < bestDist) { bestDist = dist; best = j; }
            }
            Assert.True(best >= 0 && bestDist <= tol,
                $"eigenvalue {e} has no match within {tol:E1} (nearest remaining at {bestDist:E2})");
            remaining.RemoveAt(best);
        }
        Assert.Empty(remaining);
    }
}
