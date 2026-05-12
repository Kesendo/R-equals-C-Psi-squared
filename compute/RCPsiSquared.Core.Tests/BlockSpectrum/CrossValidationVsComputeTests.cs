using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using ComputeLiouvillian = RCPsiSquared.Compute.Liouvillian;
using ComputeTopology = RCPsiSquared.Compute.Topology;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

/// <summary>Phase A cross-validation: assert that <see cref="LiouvillianBlockSpectrum"/>
/// (Core: <see cref="PauliHamiltonian.XYChain"/> + <see cref="PauliDephasingDissipator.BuildZ"/>
/// + per-block eig over joint-popcount sectors) and <see cref="ComputeLiouvillian"/>
/// (Compute: <see cref="ComputeTopology.ChainXY"/> + <see cref="ComputeLiouvillian.BuildDirectRaw"/>
/// + <see cref="ComputeLiouvillian.GetAllEigenvaluesMklRaw"/>) produce bit-exactly the same
/// spectrum as a multiset for the uniform XY chain + uniform Z-dephasing Liouvillian.
///
/// <para>Both pipelines build the same physical L:
/// <c>L = −i [H, ·] + Σ_l γ_l (Z_l ρ Z_l − ρ)</c> with
/// <c>H = (J/2) Σ_b (X_b X_{b+1} + Y_b Y_{b+1})</c> on a chain of N sites, uniform
/// <c>γ_l = γ</c>. Convention pairs verified equivalent across the two implementations:
/// <list type="bullet">
///   <item>Hamiltonian factor: both use <c>J/2</c> per (XX+YY) bilinear (Core via
///         <c>XYChain</c>; Compute via <c>ChainXY</c>'s <c>J/2</c> halving).</item>
///   <item>Dissipator factor: both use rate <c>γ_l</c> with Lindblad operator <c>√γ_l Z_l</c>
///         (Core: <c>γ (P ρ P − ρ)</c>; Compute: <c>L_k = √γ Z_k</c> ⇒ <c>γ Z ρ Z − γ ρ</c>),
///         identical generators since <c>Z² = I</c>.</item>
///   <item>Vec layout: both row-major <c>flat = i·d + j</c> (Compute writes the matching
///         column-major LAPACK layout <c>data[col·d² + row]</c> — same matrix entries).</item>
///   <item>Commutator sign: both <c>−i (H ⊗ I − I ⊗ H^T)</c>.</item>
///   <item>Qubit order: both site 0 = leftmost Kronecker factor = MSB.</item>
/// </list>
/// Eigenvalues are basis-invariant, so identical L matrices give identical multisets
/// regardless of internal storage order.</para>
///
/// <para>Anchors <see cref="LiouvillianBlockSpectrum"/> against the Compute project as
/// ground truth before extending BlockSpectrum to broader Hamiltonians, dissipators, and
/// topologies (Phase B+).</para></summary>
public class CrossValidationVsComputeTests
{
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void Spectrum_BlockSpectrum_Matches_Compute_FullEig_XYChain(int N)
    {
        const double J = 1.0;
        const double gamma = 0.5;
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();

        // Pipeline 1: Core BlockSpectrum (per-block eig over joint-popcount sectors,
        // never materialises the full 4^N × 4^N L).
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        Complex[] blockSpectrum = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N);

        // Pipeline 2: Compute ground truth (full dense L → direct LAPACK eigvals via MKL).
        var couplings = Enumerable.Repeat(J, N - 1).ToArray();
        var bonds = ComputeTopology.ChainXY(N, couplings);
        Complex[] lRaw = ComputeLiouvillian.BuildDirectRaw(N, bonds, gammaPerSite);
        int liouvilleDim = 1 << (2 * N);
        Complex[] computeSpectrum = ComputeLiouvillian.GetAllEigenvaluesMklRaw(lRaw, liouvilleDim);

        Assert.Equal(liouvilleDim, blockSpectrum.Length);
        Assert.Equal(liouvilleDim, computeSpectrum.Length);
        AssertMultisetEqual(blockSpectrum, computeSpectrum, tol: 1e-9, N: N);
    }

    /// <summary>Greedy nearest-neighbour multiset comparison, robust to degenerate
    /// eigenvalue clusters (XY chain spectra are highly degenerate under joint-popcount
    /// blocking — a sort-then-zip would mis-pair conjugate-pair members across the two
    /// independent eigensolvers). Mirrors <see cref="LiouvillianBlockSpectrumTests"/>
    /// for consistency.</summary>
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
