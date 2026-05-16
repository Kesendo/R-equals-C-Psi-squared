using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using Xunit;
using Xunit.Abstractions;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.BlockSpectrum.JordanWigner;

/// <summary>Witnesses for <see cref="JwSlaterPairArnoldiEig"/>: managed Arnoldi iteration
/// on the CSR sparse L_JW from <see cref="JwSlaterPairSparseLBuilder"/>. Validated against
/// MathNet's dense Evd at small N (4–6); demonstrates top-k eigenvalue extraction at the
/// N=10 max-block where dense Evd is infeasible.</summary>
public class JwSlaterPairArnoldiEigTests
{
    private readonly ITestOutputHelper _out;
    public JwSlaterPairArnoldiEigTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(4, 1, 2, 24, 1e-6)]   // dim ≤ numIter → Arnoldi converges to FP roundoff (~1e-7)
    [InlineData(5, 2, 3, 70, 1e-7)]   // m=70 / dim 100 ≈ 70%: tight relative convergence
    [InlineData(6, 3, 3, 80, 1e-3)]   // m=80 / dim 400 = 20%: looser practical bound
    public void Build_LargestMagnitude_MatchesDenseEvd(int N, int pCol, int pRow, int numIter, double relTolerance)
    {
        // Largest-magnitude Ritz values match the corresponding dense Evd eigenvalues to a
        // relative tolerance that loosens as the Krylov fraction (numIter / sectorDim)
        // shrinks. Plain Arnoldi without implicit restart needs Krylov dim of substantial
        // fraction of the spectrum to converge tightly across the requested numEig values.
        var gamma = Enumerable.Repeat(0.1, N).ToArray();
        var sparse = JwSlaterPairSparseLBuilder.Build(N, pCol, pRow, gamma);

        int dim = sparse.SectorDim;
        int numEig = Math.Min(8, dim - 1);
        int actualNumIter = Math.Min(numIter, dim - 1);
        var arnoldi = JwSlaterPairArnoldiEig.Build(sparse, numEig, actualNumIter, randomSeed: 42);

        var dense = sparse.ToDense();
        var denseEigs = dense.Evd().EigenValues.ToArray()
            .OrderByDescending(z => z.Magnitude).Take(numEig).ToArray();

        double maxMatchDistance = 0.0;
        foreach (var arnEig in arnoldi.Eigenvalues)
        {
            double nearest = denseEigs.Min(d => (d - arnEig).Magnitude);
            if (nearest > maxMatchDistance) maxMatchDistance = nearest;
        }
        double scale = denseEigs.Max(d => d.Magnitude);
        _out.WriteLine($"N={N}, (p_c,p_r)=({pCol},{pRow}), dim={dim}, numEig={numEig}, " +
                       $"numIter={actualNumIter}: max-match-distance={maxMatchDistance:G3}, " +
                       $"scale={scale:F3}, relative={maxMatchDistance / scale:G3}, tol={relTolerance:G3}");
        Assert.True(maxMatchDistance / scale < relTolerance,
            $"Max relative match distance {maxMatchDistance / scale:G3} exceeds tolerance {relTolerance:G3}");
    }

    [Theory]
    [InlineData(4, 2, 2, 30, 1e-7)]   // half-filling: m=30 / dim 36 ≈ 83%, tight convergence
    [InlineData(6, 3, 3, 200, 1e-3)]  // half-filling at N=6: m=200 / dim 400 = 50%, practical bound
    public void Build_PalindromicGamma_RecoversPalindromeCenter(int N, int pCol, int pRow, int numIter, double mirrorTol)
    {
        // F1 palindrome holds for the FULL Liouvillian spectrum, paired around −Σγ. Per-
        // sector it holds only at half-filling (p_c = p_r = N/2) — other sectors are paired
        // with their X⊗N partner (N−p_c, N−p_r). For each Arnoldi-recovered eigenvalue λ,
        // the mirror λ' = −2·Σγ − λ must appear in the full sector spectrum.
        // Mirror-distance tolerance tracks Arnoldi convergence at the chosen Krylov fraction.
        if (pCol != N / 2 || pRow != N / 2)
            throw new InvalidOperationException("This test asserts F1 palindrome at half-filling only.");

        double gamma = 0.07;
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();
        double center = -N * gamma;

        var sparse = JwSlaterPairSparseLBuilder.Build(N, pCol, pRow, gammaPerSite);
        var arnoldi = JwSlaterPairArnoldiEig.Build(sparse, numEig: 6, numIter: numIter, randomSeed: 7);

        var dense = sparse.ToDense();
        var allEigs = dense.Evd().EigenValues.ToArray();

        double maxMirrorDistance = 0.0;
        foreach (var arnEig in arnoldi.Eigenvalues)
        {
            var mirror = new Complex(2 * center - arnEig.Real, -arnEig.Imaginary);
            double nearest = allEigs.Min(e => (e - mirror).Magnitude);
            if (nearest > maxMirrorDistance) maxMirrorDistance = nearest;
        }
        _out.WriteLine($"N={N}, half-filling (p_c=p_r={pCol}), center=−Σγ={center}, numIter={numIter}: " +
                       $"max-mirror-distance={maxMirrorDistance:G3}, tol={mirrorTol:G3}");
        Assert.True(maxMirrorDistance < mirrorTol,
            $"Palindrome mirror distance {maxMirrorDistance:G3} exceeds tolerance {mirrorTol:G3} — F1 violated or Arnoldi under-converged?");
    }

    [Fact]
    public void Build_Tier_IsTier1Derived()
    {
        var gamma = Enumerable.Repeat(0.1, 4).ToArray();
        var sparse = JwSlaterPairSparseLBuilder.Build(N: 4, pCol: 1, pRow: 2, gamma);
        var arnoldi = JwSlaterPairArnoldiEig.Build(sparse, numEig: 4, numIter: 16, randomSeed: 1);
        Assert.Equal(Tier.Tier1Derived, arnoldi.Tier);
    }

    [Fact]
    public void Build_RejectsNumEigExceedingNumIter()
    {
        var gamma = Enumerable.Repeat(0.1, 4).ToArray();
        var sparse = JwSlaterPairSparseLBuilder.Build(N: 4, pCol: 1, pRow: 2, gamma);
        Assert.Throws<ArgumentException>(() =>
            JwSlaterPairArnoldiEig.Build(sparse, numEig: 10, numIter: 5, randomSeed: 1));
    }

    [Fact]
    public void Build_AtN10HalfFilling_ExtractsLargestMagnitudeEigenvalues()
    {
        // N=10 (5,5) max-block: dim 63504, no dense Evd possible. Compute top-8 by Arnoldi
        // and verify (a) they're returned without OOM, (b) magnitudes are bounded by the
        // operator's Frobenius norm. The witness is "we can extract eigenvalues at N=10".
        // γ₀ = 0.05 = F86 sweep convention; the absolute |λ| values reported are γ-dependent
        // sample artifacts, not framework-canonical numbers (PRIMORDIAL_GAMMA_CONSTANT:
        // γ₀ universal but unspecified; only Q = J/γ₀ is intrinsic).
        double gamma = 0.05;
        var gammaPerSite = Enumerable.Repeat(gamma, 10).ToArray();
        var sparse = JwSlaterPairSparseLBuilder.Build(N: 10, pCol: 5, pRow: 5, gammaPerSite);

        var arnoldi = JwSlaterPairArnoldiEig.Build(sparse, numEig: 8, numIter: 40, randomSeed: 1);

        // Frobenius² of L_JW = Σ|values|² (CSR). Each eigenvalue must satisfy |λ| ≤ ‖L‖_F.
        double frobSq = sparse.Values.Sum(v => v.Real * v.Real + v.Imaginary * v.Imaginary);
        double frobNorm = Math.Sqrt(frobSq);

        foreach (var e in arnoldi.Eigenvalues)
            Assert.True(e.Magnitude <= frobNorm,
                $"Eigenvalue |λ|={e.Magnitude:F3} exceeds ‖L‖_F={frobNorm:F3}");

        _out.WriteLine($"N=10 (5,5) sectorDim={sparse.SectorDim}: ‖L‖_F={frobNorm:F3}, " +
                       $"top-8 |λ|: " +
                       string.Join(", ", arnoldi.Eigenvalues.Select(e => e.Magnitude.ToString("F3"))));
    }
}
