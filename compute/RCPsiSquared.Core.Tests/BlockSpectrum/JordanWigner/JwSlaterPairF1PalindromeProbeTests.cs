using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.BlockSpectrum.JordanWigner;

/// <summary>Witnesses for <see cref="JwSlaterPairF1PalindromeProbe"/>: two shift-invert
/// Arnoldi extractions at σ_slow ≈ 0 and σ_fast = −2·Σγ − σ_slow, per-pair residuals on the
/// F1 mirror map <c>λ → −λ − 2·Σγ</c> with the L-conjugation flexibility absorbed.
///
/// <para>The probe is informational, not a strict F1 prover (see the class docstring).
/// Tests therefore assert structural properties — extraction completes, returned eigenvalues
/// genuinely lie in the spectrum, σ_fast equals the F1 mirror of σ_slow exactly — and the
/// numerical "match count" at a coarse tolerance for the K-pair sample. Tight per-pair
/// matching depends on the slow and fast Arnoldi runs landing on F1-mirrored Krylov
/// subspaces, which is governed by Arnoldi depth × spectrum density and is not guaranteed
/// at finite numIter.</para></summary>
public class JwSlaterPairF1PalindromeProbeTests
{
    private readonly ITestOutputHelper _out;
    public JwSlaterPairF1PalindromeProbeTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(5, 2, 2)]
    [InlineData(6, 3, 3)]
    [InlineData(6, 2, 3)]
    public void Build_SmallN_RecoveredEigenvaluesLieInDenseSpectrum(int N, int pCol, int pRow)
    {
        // The fundamental sanity check on the probe: every eigenvalue returned by either end's
        // shift-invert Arnoldi must be a genuine eigenvalue of the sparse block (cross-checked
        // against dense Evd). This is independent of whether the slow/fast pairs F1-match
        // each other — it certifies the probe doesn't return Arnoldi hallucinations.
        var gamma = Enumerable.Repeat(0.05, N).ToArray();
        var sparse = JwSlaterPairSparseLBuilder.Build(N, pCol, pRow, gamma);
        int dim = sparse.SectorDim;
        int numEig = 4;
        int numIter = Math.Min(30, dim - 1);

        var probe = JwSlaterPairF1PalindromeProbe.Build(sparse,
            numEig, numIter, randomSeed: 17, innerTolerance: 1e-10, innerMaxIter: 500);

        var dense = sparse.ToDense();
        var denseEigs = dense.Evd().EigenValues.ToArray();
        foreach (var slow in probe.SlowExtraction.Eigenvalues)
        {
            double dMin = denseEigs.Min(e => (e - slow).Magnitude);
            Assert.True(dMin < 1e-6,
                $"Slow eigenvalue ({slow.Real}, {slow.Imaginary}) not in dense spectrum (min dist {dMin:G3})");
        }
        foreach (var fast in probe.FastExtraction.Eigenvalues)
        {
            double dMin = denseEigs.Min(e => (e - fast).Magnitude);
            Assert.True(dMin < 1e-6,
                $"Fast eigenvalue ({fast.Real}, {fast.Imaginary}) not in dense spectrum (min dist {dMin:G3})");
        }

        _out.WriteLine($"N={N} ({pCol},{pRow}), dim={dim}, 2·Σγ={probe.SumGamma * 2:F4}: " +
                       $"max F1 residual = {probe.MaxF1Residual:G3}, mean = {probe.MeanF1Residual:G3}, " +
                       $"matched(@1e-4)={probe.CountMatchedPairs(1e-4)}/{numEig}");
    }

    [Fact]
    public void Build_Tier_IsTier1Derived()
    {
        var gamma = Enumerable.Repeat(0.1, 4).ToArray();
        var sparse = JwSlaterPairSparseLBuilder.Build(N: 4, pCol: 2, pRow: 2, gamma);
        var probe = JwSlaterPairF1PalindromeProbe.Build(sparse,
            numEig: 3, numIter: 15, randomSeed: 1, innerTolerance: 1e-10, innerMaxIter: 200);
        Assert.Equal(Tier.Tier1Derived, probe.Tier);
    }

    [Fact]
    public void Build_MirrorShiftRelationshipIsExact()
    {
        // σ_fast must equal −σ_slow − 2·Σγ exactly (the F1 mirror of σ_slow), so the two
        // shift-invert runs operate at matched conditioning. This is the structural property
        // the probe contracts to its caller, independent of Arnoldi convergence quality.
        var gamma = Enumerable.Repeat(0.1, 4).ToArray();
        var sparse = JwSlaterPairSparseLBuilder.Build(N: 4, pCol: 2, pRow: 2, gamma);

        var sigmaSlow = new Complex(0.01, 0.005);
        var probe = JwSlaterPairF1PalindromeProbe.Build(sparse,
            numEig: 3, numIter: 15, randomSeed: 1, innerTolerance: 1e-10, innerMaxIter: 200,
            sigmaSlow: sigmaSlow);

        double sumGamma = gamma.Sum();
        var expectedFast = new Complex(-2 * sumGamma - sigmaSlow.Real, -sigmaSlow.Imaginary);
        Assert.Equal(expectedFast.Real, probe.SigmaFast.Real, 12);
        Assert.Equal(expectedFast.Imaginary, probe.SigmaFast.Imaginary, 12);
    }

    [Fact]
    public void Build_AtN10HalfFilling_RunsAndRecoveredEigenvaluesAreBounded()
    {
        // The N=10 (5,5) max-block: dim 63504, dense Evd infeasible (~64 GB). This is the
        // probe's target case — it cannot be cross-validated against dense Evd at this size.
        // The structural assertions: probe completes within reasonable time, both extractions
        // return eigenvalues bounded by ‖L‖_F (no Arnoldi blowup), and a brief reconnaissance
        // log of the slow/fast pairs + per-pair residuals is emitted for human inspection.
        var gamma = Enumerable.Repeat(0.05, 10).ToArray();

        var sw = System.Diagnostics.Stopwatch.StartNew();
        var sparse = JwSlaterPairSparseLBuilder.Build(N: 10, pCol: 5, pRow: 5, gamma);
        sw.Stop();
        var sparseBuildElapsed = sw.Elapsed;
        _out.WriteLine($"sparse L_JW build at N=10 (5,5): dim={sparse.SectorDim}, nnz={sparse.NnzTotal:N0}, " +
                       $"max-nnz/row={sparse.MaxNnzPerRow}, build = {sparseBuildElapsed.TotalSeconds:F1} s");

        sw.Restart();
        var probe = JwSlaterPairF1PalindromeProbe.Build(sparse,
            numEig: 4, numIter: 20, randomSeed: 1, innerTolerance: 1e-8, innerMaxIter: 1000);
        sw.Stop();

        _out.WriteLine($"two-shift probe elapsed = {sw.Elapsed.TotalSeconds:F1} s, " +
                       $"2·Σγ = {probe.SumGamma * 2:F4}");
        _out.WriteLine($"σ_slow = {probe.SigmaSlow}, σ_fast = {probe.SigmaFast}");
        _out.WriteLine($"slow inner BiCGStab mean/max iter = " +
                       $"{probe.SlowExtraction.MeanInnerIterations:F1} / {probe.SlowExtraction.MaxInnerIterations}");
        _out.WriteLine($"fast inner BiCGStab mean/max iter = " +
                       $"{probe.FastExtraction.MeanInnerIterations:F1} / {probe.FastExtraction.MaxInnerIterations}");

        for (int i = 0; i < probe.SlowExtraction.Eigenvalues.Length; i++)
        {
            var s = probe.SlowExtraction.Eigenvalues[i];
            var pred = new Complex(-s.Real - 2 * probe.SumGamma, -s.Imaginary);
            _out.WriteLine($"  slow_{i} = ({s.Real:F6}, {s.Imaginary:F6}) → predicted mirror " +
                           $"({pred.Real:F6}, {pred.Imaginary:F6}), residual = {probe.PerSlowResidual[i]:G3}");
        }
        foreach (var f in probe.FastExtraction.Eigenvalues)
            _out.WriteLine($"  fast    = ({f.Real:F6}, {f.Imaginary:F6})");
        _out.WriteLine($"max F1 residual = {probe.MaxF1Residual:G3}, mean = {probe.MeanF1Residual:G3}, " +
                       $"matched(@1e-4)={probe.CountMatchedPairs(1e-4)}/4");

        // Structural assertion: all returned eigenvalues bounded by sparse-Frobenius-norm
        // (no Arnoldi blow-up; no NaN).
        double frobNorm = Math.Sqrt(sparse.Values.Sum(v => v.Real * v.Real + v.Imaginary * v.Imaginary));
        foreach (var e in probe.SlowExtraction.Eigenvalues.Concat(probe.FastExtraction.Eigenvalues))
            Assert.True(e.Magnitude <= frobNorm,
                $"Eigenvalue |λ|={e.Magnitude:F3} exceeds ‖L‖_F={frobNorm:F3}");
    }
}
