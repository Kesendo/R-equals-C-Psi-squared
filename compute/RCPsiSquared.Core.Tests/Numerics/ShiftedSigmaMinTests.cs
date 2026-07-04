using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.Numerics;

/// <summary>The step-3 census's σ_min probe: one dense LU (zgetrf) + inverse power iteration on
/// [(M−s)ᴴ(M−s)]⁻¹, no full spectrum. Cross-validated against MathNet dense SVD on a random complex
/// matrix and on a physical coherence block at the N=5 seed locus; the exactly-singular-shift path
/// must degrade to σ_min = 0 without crashing.</summary>
public class ShiftedSigmaMinTests
{
    static double DenseSigmaMin(Complex[,] m, Complex shift)
    {
        int d = m.GetLength(0);
        var dense = Matrix<Complex>.Build.DenseOfArray(m)
                    - Matrix<Complex>.Build.DenseIdentity(d) * shift;
        return dense.Svd().S.Enumerate().Select(z => z.Magnitude).Min();
    }

    [Fact]
    public void MatchesDenseSvd_OnRandomComplexMatrix()
    {
        var rng = new Random(7);
        int d = 60;
        var m = new Complex[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                m[i, j] = new Complex(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
        var shift = new Complex(0.3, -0.2);

        var probe = ShiftedSigmaMin.Estimate(m, shift, maxIter: 500);
        double svd = DenseSigmaMin(m, shift);

        Assert.True(probe.Converged);
        Assert.Equal(svd, probe.SigmaMin, svd * 1e-4);
    }

    [Fact]
    public void MatchesDenseSvd_OnPhysicalBlock()
    {
        // the (2,3) block at the N=5 seed locus; shift = the recorded lambda_A
        var m = WeightCoherenceBlock.Build(5, 2, 3, new Complex(0.620878, 0));
        var shift = new Complex(-4.618886, 0);
        var probe = ShiftedSigmaMin.Estimate(m, shift, maxIter: 500);
        double svd = DenseSigmaMin(m, shift);
        Assert.True(probe.Converged);
        Assert.Equal(svd, probe.SigmaMin, Math.Max(svd * 1e-3, 1e-10));
    }

    [Fact]
    public void ExactEigenvalueShift_ReportsTinySigma_WithoutCrashing()
    {
        // shift = an exact computed eigenvalue of a small block
        var m = WeightCoherenceBlock.Build(4, 1, 1, new Complex(0.7, 0));
        var ev = Matrix<Complex>.Build.DenseOfArray(m).Evd().EigenValues[0];
        var probe = ShiftedSigmaMin.Estimate(m, ev);
        Assert.True(probe.SigmaMin < 1e-10);
    }
}
