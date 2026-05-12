using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Pauli;
using Xunit;
using ComputeLiouvillian = RCPsiSquared.Compute.Liouvillian;
using ComputeTopology = RCPsiSquared.Compute.Topology;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.SymmetryFamily;

/// <summary>Phase 7 final cross-validation: anchors F92 (anti-palindromic J spectral
/// invariance) against the Compute project's full-LAPACK ground-truth pipeline. Builds
/// the same physical Liouvillian L two ways:
///
/// <list type="number">
///   <item><b>Pipeline 1 (Core BlockSpectrum):</b>
///         <see cref="PauliHamiltonian"/> from explicit per-bond (J/2)·(X_b X_{b+1} + Y_b Y_{b+1})
///         <c>PauliTerm</c>s + <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/>
///         (per-block eig over joint-popcount sectors, full L spectrum).</item>
///   <item><b>Pipeline 2 (Compute ground truth):</b>
///         <see cref="ComputeTopology.ChainXY"/> with per-bond <c>couplings</c> array +
///         <see cref="ComputeLiouvillian.BuildDirectRaw"/> + <see cref="ComputeLiouvillian.GetAllEigenvaluesMklRaw"/>
///         (full dense L → direct LAPACK).</item>
/// </list>
///
/// <para>Asserts bit-exact multiset match (|Δλ| &lt; 1e-9) at N=4, 5 for the same
/// anti-palindromic linear-gradient J profile used by
/// <c>F92BondAntiPalindromicJSpectralInvarianceTests</c>. This is the FULL spectrum
/// (not the F71-refined diagonal block): the F92 invariance test verifies that the
/// F71-refined diagonal-block spectrum is invariant under anti-palindromic J; THIS test
/// verifies that two independent pipelines agree on the full L spectrum for the same
/// inhomogeneous J. Together they nail down: (a) BlockSpectrum reproduces ground truth
/// even with inhomogeneous J, (b) F71-refined sub-spectrum is invariant under
/// anti-palindromic perturbations of that J. Closes the SymmetryFamily Phase 7
/// implementation arc.</para></summary>
public sealed class CrossValidationVsComputeTests
{
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    public void F92_AntiPalindromicJ_BlockSpectrum_Matches_Compute_FullEig(int N)
    {
        const double gamma = 0.5;
        const double Javg = 1.0;
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();

        // Build anti-palindromic J profile: J_b + J_{N-2-b} = 2·Javg via linear gradient
        // J_b = Javg + 0.3·(b - (N-2)/2). Same profile as F92 invariance test.
        int numBonds = N - 1;
        var bondJ = new double[numBonds];
        double mid = (numBonds - 1) / 2.0;
        for (int b = 0; b < numBonds; b++) bondJ[b] = Javg + 0.3 * (b - mid);

        // Pipeline 1: Core BlockSpectrum (per-block eig over joint-popcount sectors,
        // full L spectrum without ever materialising 4^N × 4^N L).
        var allTerms = new List<PauliTerm>();
        for (int b = 0; b < numBonds; b++)
        {
            Complex c = bondJ[b] / 2.0;  // F86 convention (J/2) per PauliHamiltonian.XYChain
            allTerms.Add(PauliTerm.TwoSite(N, b, PauliLetter.X, b + 1, PauliLetter.X, c));
            allTerms.Add(PauliTerm.TwoSite(N, b, PauliLetter.Y, b + 1, PauliLetter.Y, c));
        }
        ComplexMatrix H = new PauliHamiltonian(N, allTerms).ToMatrix();
        Complex[] blockSpectrum = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N);

        // Pipeline 2: Compute ground truth (full dense L → direct LAPACK eigvals via MKL).
        var bonds = ComputeTopology.ChainXY(N, bondJ);
        Complex[] lRaw = ComputeLiouvillian.BuildDirectRaw(N, bonds, gammaPerSite);
        int liouvilleDim = 1 << (2 * N);
        Complex[] computeSpectrum = ComputeLiouvillian.GetAllEigenvaluesMklRaw(lRaw, liouvilleDim);

        Assert.Equal(liouvilleDim, blockSpectrum.Length);
        Assert.Equal(liouvilleDim, computeSpectrum.Length);
        AssertMultisetEqual(blockSpectrum, computeSpectrum, tol: 1e-9, N: N);
    }

    /// <summary>Greedy nearest-neighbour multiset comparison, robust to degenerate
    /// eigenvalue clusters (XY chain spectra are highly degenerate under joint-popcount
    /// blocking; sort-then-zip would mis-pair conjugate-pair members across the two
    /// independent eigensolvers). Mirrors the Phase A pattern in
    /// <c>BlockSpectrum.CrossValidationVsComputeTests</c>.</summary>
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
                $"N={N}: eigenvalue {i} mismatch, expected={expected[i]}, nearest actual={(bestJ >= 0 ? actual[bestJ].ToString() : "<none>")}, |Δ|={bestDist:E3} (tol={tol:E1}).");
            taken[bestJ] = true;
        }
    }
}
