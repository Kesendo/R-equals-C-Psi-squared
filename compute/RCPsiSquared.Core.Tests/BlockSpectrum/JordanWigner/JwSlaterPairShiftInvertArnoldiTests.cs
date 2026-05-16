using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using Xunit;
using Xunit.Abstractions;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.BlockSpectrum.JordanWigner;

/// <summary>Witnesses for <see cref="JwSlaterPairShiftInvertArnoldi"/>: shift-invert
/// Arnoldi on the CSR sparse L_JW. Computes eigenvalues nearest a complex shift σ by
/// running plain Arnoldi on the operator (L − σI)^(−1); the inner linear solve uses an
/// in-house BiCGStab on the same CSR matvec. Validated against dense Evd at small N;
/// extracts slow modes (Re λ ≈ 0) at the N=10 max-block where direct dense eig is
/// infeasible.</summary>
public class JwSlaterPairShiftInvertArnoldiTests
{
    private readonly ITestOutputHelper _out;
    public JwSlaterPairShiftInvertArnoldiTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(4, 1, 2, 0.0, 1e-7)]   // dim 24, m=20 ≈ 83% Krylov: tight
    [InlineData(4, 2, 2, 0.0, 1e-7)]   // dim 36, m=20 ≈ 55%: tight
    [InlineData(5, 2, 3, 0.0, 1e-2)]   // dim 100, m=20 = 20%: plain shift-invert practical bound
    public void Build_SigmaNearZero_RecoversSlowModes(int N, int pCol, int pRow, double sigmaReal, double matchTol)
    {
        // For shift σ slightly off the steady state, shift-invert Arnoldi extracts the
        // eigenvalues nearest σ. With σ ≈ 0 we get the slow modes (Re λ ≈ 0) — the
        // physically interesting decoherence-time-scale modes.
        // Sigma offset: use σ = sigmaReal + 0.001i to dodge the L = 0 steady state.
        var sigma = new Complex(sigmaReal, 0.001);
        var gamma = Enumerable.Repeat(0.1, N).ToArray();
        var sparse = JwSlaterPairSparseLBuilder.Build(N, pCol, pRow, gamma);

        int dim = sparse.SectorDim;
        int numEig = Math.Min(4, dim - 2);
        int numIter = Math.Min(20, dim - 1);

        var shiftInvert = JwSlaterPairShiftInvertArnoldi.Build(sparse, sigma,
            numEig, numIter, randomSeed: 13, innerTolerance: 1e-10, innerMaxIter: 500);

        var dense = sparse.ToDense();
        var allEigs = dense.Evd().EigenValues.ToArray();
        var nearestSigma = allEigs.OrderBy(e => (e - sigma).Magnitude).Take(numEig).ToArray();

        double maxMatch = 0.0;
        foreach (var ev in shiftInvert.Eigenvalues)
        {
            double nearest = nearestSigma.Min(d => (d - ev).Magnitude);
            if (nearest > maxMatch) maxMatch = nearest;
        }
        _out.WriteLine($"N={N}, ({pCol},{pRow}), dim={dim}, σ={sigma}, numEig={numEig}, " +
                       $"numIter={numIter}: max-match-to-nearest-σ = {maxMatch:G3}, tol = {matchTol:G3}");
        _out.WriteLine($"  Recovered: {string.Join(", ", shiftInvert.Eigenvalues.Select(e => $"({e.Real:F4},{e.Imaginary:F4})"))}");
        _out.WriteLine($"  Expected:  {string.Join(", ", nearestSigma.Select(e => $"({e.Real:F4},{e.Imaginary:F4})"))}");
        Assert.True(maxMatch < matchTol,
            $"Max distance {maxMatch:G3} from nearest-σ exceeds tolerance {matchTol:G3}");
    }

    [Fact]
    public void Build_Tier_IsTier1Derived()
    {
        var gamma = Enumerable.Repeat(0.1, 4).ToArray();
        var sparse = JwSlaterPairSparseLBuilder.Build(N: 4, pCol: 1, pRow: 2, gamma);
        var shiftInvert = JwSlaterPairShiftInvertArnoldi.Build(sparse, new Complex(0, 0.001),
            numEig: 4, numIter: 10, randomSeed: 1, innerTolerance: 1e-10, innerMaxIter: 200);
        Assert.Equal(Tier.Tier1Derived, shiftInvert.Tier);
    }

    [Fact]
    public void Build_RejectsNumEigExceedingNumIter()
    {
        var gamma = Enumerable.Repeat(0.1, 4).ToArray();
        var sparse = JwSlaterPairSparseLBuilder.Build(N: 4, pCol: 1, pRow: 2, gamma);
        Assert.Throws<ArgumentException>(() =>
            JwSlaterPairShiftInvertArnoldi.Build(sparse, new Complex(0, 0.001),
                numEig: 10, numIter: 5, randomSeed: 1, innerTolerance: 1e-10, innerMaxIter: 100));
    }

    [Fact]
    public void Build_AtN10HalfFilling_ExtractsSlowModesNearZero()
    {
        // N=10 (5,5) max-block, σ near 0: extract the 6 slowest modes (Re λ closest to 0).
        // γ₀ = 0.05 = F86 sweep convention; absolute eigenvalue magnitudes reported below
        // are γ-dependent sample artifacts (PRIMORDIAL_GAMMA_CONSTANT: γ₀ universal but
        // unspecified; only Q = J/γ₀ is intrinsic).
        // No dense reference available; the witness is "Arnoldi converges and the extracted
        // eigenvalues are sub-magnitude of the operator Frobenius norm".
        var gamma = Enumerable.Repeat(0.05, 10).ToArray();
        var sparse = JwSlaterPairSparseLBuilder.Build(N: 10, pCol: 5, pRow: 5, gamma);

        var sigma = new Complex(0, 0.001);
        var shiftInvert = JwSlaterPairShiftInvertArnoldi.Build(sparse, sigma,
            numEig: 6, numIter: 20, randomSeed: 1, innerTolerance: 1e-8, innerMaxIter: 1000);

        double frobNorm = Math.Sqrt(sparse.Values.Sum(v => v.Real * v.Real + v.Imaginary * v.Imaginary));
        foreach (var e in shiftInvert.Eigenvalues)
            Assert.True(e.Magnitude <= frobNorm,
                $"Slow eigenvalue |λ|={e.Magnitude:F3} exceeds ‖L‖_F={frobNorm:F3}");

        // Each "slow" eigenvalue should have Re λ close to σ_real = 0; bound generously by
        // Σγ (= 0.5 at γ₀=0.05, N=10) — slow modes near 0 sit well within this band.
        double sumGamma = gamma.Sum();
        foreach (var e in shiftInvert.Eigenvalues)
            Assert.True(Math.Abs(e.Real) < sumGamma,
                $"Slow eigenvalue Re={e.Real:F4} should be near 0 (σ_real); shift-invert may have under-converged");

        _out.WriteLine($"N=10 (5,5) γ₀=0.05, σ={sigma}: slow eigenvalues (sample artifact, scales with γ₀):");
        foreach (var e in shiftInvert.Eigenvalues)
            _out.WriteLine($"  λ = ({e.Real:F6}, {e.Imaginary:F6}),  |λ| = {e.Magnitude:F6}");
        _out.WriteLine($"  inner-solver mean iter = {shiftInvert.MeanInnerIterations:F1}, " +
                       $"max = {shiftInvert.MaxInnerIterations}");
    }
}
