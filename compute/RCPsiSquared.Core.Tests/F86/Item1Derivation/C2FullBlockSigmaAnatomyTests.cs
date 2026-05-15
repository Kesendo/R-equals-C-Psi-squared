using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class C2FullBlockSigmaAnatomyTests : IClassFixture<C2FullBlockSigmaAnatomyCache>
{
    private readonly ITestOutputHelper _out;
    private readonly C2FullBlockSigmaAnatomyCache _cache;
    public C2FullBlockSigmaAnatomyTests(C2FullBlockSigmaAnatomyCache cache, ITestOutputHelper @out)
    {
        _cache = cache;
        _out = @out;
    }

    private static CoherenceBlock C2Block(int N) =>
        new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);

    [Fact]
    public void Build_AtC2N5_ReturnsTier2VerifiedClaim()
    {
        var anatomy = _cache.Get(5);
        Assert.Equal(Tier.Tier2Verified, anatomy.Tier);
        Assert.Equal(5, anatomy.Block.N);
    }

    [Fact]
    public void Build_AtC2N5_ProducesOneWitnessPerEigenmode()
    {
        var anatomy = _cache.Get(5);
        Assert.Equal(anatomy.Block.Basis.MTotal, anatomy.SigmaSpectrum.Count);
    }

    [Fact]
    public void SigmaWitness_AtC2N5_HasNonNegativeSigma()
    {
        var anatomy = _cache.Get(5);
        Assert.All(anatomy.SigmaSpectrum, w =>
            Assert.True(w.Sigma >= -1e-12,
                $"Sigma must be non-negative (S is PSD); got {w.Sigma} at λ={w.EigenvalueReal}+{w.EigenvalueImag}i"));
    }

    [Fact]
    public void SigmaWitness_AtC2N5_HasPositiveTotalSigma()
    {
        var anatomy = _cache.Get(5);
        double total = anatomy.SigmaSpectrum.Sum(w => w.Sigma);
        Assert.True(total > 1e-6, $"Total sigma must be positive; got {total}");
    }

    [Theory]
    [InlineData(4, 2)]    // c=2 N=4: path-3, smallest case, no bare site
    [InlineData(5, 2)]    // c=2 N=5: path-4 → F_a count = floor(5/2) = 2
    [InlineData(6, 3)]    // c=2 N=6: path-5 → F_a count = floor(6/2) = 3
    [InlineData(7, 3)]    // c=2 N=7: path-6 → F_a count = floor(7/2) = 3
    [InlineData(8, 4)]    // c=2 N=8: path-7 → F_a count = floor(8/2) = 4
    [InlineData(10, 5)]   // c=2 N=10: path-9, beyond F89 closed-form table; stress check
    public void FaModes_Count_MatchesFaCount(int n, int expectedFaCount)
    {
        var anatomy = _cache.Get(n);
        int actualFaCount = anatomy.SigmaSpectrum.Count(w => w.BlochIndexN.HasValue);
        Assert.Equal(expectedFaCount, actualFaCount);
    }


    [Fact]
    public void FaModes_BlochIndices_AreInSeAntiOrbit()
    {
        var anatomy = _cache.Get(7);   // N=7 → orbit {2, 4, 6}
        var assigned = anatomy.SigmaSpectrum
            .Where(w => w.BlochIndexN.HasValue)
            .Select(w => w.BlochIndexN!.Value)
            .OrderBy(n => n)
            .ToArray();
        Assert.Equal(new[] { 2, 4, 6 }, assigned);
    }

    [Theory]
    [InlineData(3, 2)]   // path-3 N=5: P(y_2) = 14·y_2 + 47, y_2 = 4cos(2π/5) ≈ 1.236
    [InlineData(3, 4)]
    [InlineData(4, 2)]   // path-4 N=6
    [InlineData(4, 4)]
    [InlineData(5, 2)]   // path-5 N=7
    [InlineData(5, 4)]
    [InlineData(5, 6)]
    [InlineData(6, 2)]   // path-6 N=8: includes the y_4 = 0 zero-mode case
    [InlineData(6, 4)]
    [InlineData(6, 6)]
    [InlineData(7, 2)]   // path-7 N=8: P_7(y) = 21y³ + 130y² + 292y + 382, D = 98
    [InlineData(7, 4)]
    [InlineData(7, 6)]
    [InlineData(7, 8)]
    [InlineData(8, 2)]   // path-8 N=9: P_8 = 13y³ + 54y² + 68y + 110, D = 32
    [InlineData(8, 4)]
    [InlineData(8, 6)]
    [InlineData(8, 8)]
    [InlineData(9, 2)]   // path-9 N=10: P_9 = 31y⁴ + 190y³ + 288y² + 440y + 1476, D = 324
    [InlineData(9, 4)]
    [InlineData(9, 6)]
    [InlineData(9, 8)]
    [InlineData(9, 10)]
    public void Sigma_AtPathK_MatchesF89UnifiedClosedForm(int k, int n)
    {
        int N = k + 1;   // C2Block(N) is F89 path-(N-1); for path-k take N = k+1
        var anatomy = _cache.Get(N);
        double? extracted = anatomy.SigmaForBlochIndex(n);
        Assert.NotNull(extracted);

        // Inline F89 σ formula bypassing the F89.Sigma blochN >= k+2 check
        // (we are at the boundary nBlock = N = k+1, no bare site). The math
        // is unchanged; only the runtime guard disallows this case.
        var (coefs, denom) = F89UnifiedFaClosedFormClaim.PathPolynomial(k);
        double y = F89PathKAtLockMechanismClaim.BlochEigenvalueY(k + 1, n);
        double poly = 0.0;
        double yPow = 1.0;
        foreach (var c in coefs) { poly += c * yPow; yPow *= y; }
        double expected = poly / (denom * N * N * (N - 1));

        Assert.True(Math.Abs(extracted!.Value - expected) <= 1e-8,
            $"path-{k} n={n}: extracted={extracted.Value:G10}, expected={expected:G10}, " +
            $"diff={Math.Abs(extracted.Value - expected):G6}");
    }

    [Fact]
    public void FbModes_AreInvisibleToSpatialSumKernel()
    {
        var anatomy = _cache.Get(6);   // c=2 N=6, F_b at Re(λ) ≈ -6γ₀ = -0.3
        double fbRateTarget = -6.0 * anatomy.Block.GammaZero;
        var fbModes = anatomy.SigmaSpectrum
            .Where(w => Math.Abs(w.EigenvalueReal - fbRateTarget) <= 1e-3)
            .ToList();
        Assert.NotEmpty(fbModes);
        Assert.All(fbModes, w =>
            Assert.True(w.Sigma < 1e-12,
                $"F_b mode at λ={w.EigenvalueReal}+{w.EigenvalueImag}i should have σ≈0; got {w.Sigma:G6}"));
    }

    [Fact]
    public void Sigma_AtPath7_ExtractsFourFaModes()
    {
        var anatomy = _cache.Get(8);
        var faWitnesses = anatomy.SigmaSpectrum
            .Where(w => w.BlochIndexN.HasValue)
            .OrderBy(w => w.BlochIndexN!.Value)
            .ToList();
        Assert.Equal(4, faWitnesses.Count);
        Assert.Equal(new[] { 2, 4, 6, 8 }, faWitnesses.Select(w => w.BlochIndexN!.Value).ToArray());
        Assert.All(faWitnesses, w =>
            Assert.True(w.Sigma > 0, $"path-7 F_a mode at n={w.BlochIndexN!.Value} must have σ > 0; got {w.Sigma}"));
    }

    [Fact]
    public void Sigma_AtPath7_PolynomialFitIsCubic()
    {
        var anatomy = _cache.Get(8);
        var faWitnesses = anatomy.SigmaSpectrum
            .Where(w => w.BlochIndexN.HasValue)
            .OrderBy(w => w.BlochIndexN!.Value)
            .ToList();
        Assert.Equal(4, faWitnesses.Count);

        // y_n · J for the four orbit points; σ_n · N² · (N-1) = P_path-7(y_n) (mod denom)
        int N = 8;
        int nBlock = 8;
        double[] yValues = new double[4];
        double[] sigmaScaled = new double[4];
        for (int i = 0; i < 4; i++)
        {
            var w = faWitnesses[i];
            yValues[i] = F89PathKAtLockMechanismClaim.BlochEigenvalueY(nBlock, w.BlochIndexN!.Value);
            sigmaScaled[i] = w.Sigma * N * N * (N - 1);
        }

        // Solve Vandermonde for cubic coefs c0 + c1·y + c2·y² + c3·y³ = sigmaScaled
        double[,] V = new double[4, 4];
        for (int i = 0; i < 4; i++)
        {
            double yi = yValues[i];
            V[i, 0] = 1; V[i, 1] = yi; V[i, 2] = yi * yi; V[i, 3] = yi * yi * yi;
        }
        var matV = Matrix<double>.Build.DenseOfArray(V);
        var vecB = Vector<double>.Build.DenseOfArray(sigmaScaled);
        var coefs = matV.Solve(vecB);

        // Reconstruct: P(y_i) should equal sigmaScaled[i] for each i (residual ≈ 0)
        for (int i = 0; i < 4; i++)
        {
            double yi = yValues[i];
            double reconstructed = coefs[0] + coefs[1] * yi + coefs[2] * yi * yi + coefs[3] * yi * yi * yi;
            Assert.True(Math.Abs(reconstructed - sigmaScaled[i]) <= 1e-6,
                $"Cubic fit residual at y={yi}: |reconstructed - actual| = {Math.Abs(reconstructed - sigmaScaled[i]):G6}");
        }

        _out.WriteLine($"Path-7 cubic coefficients (low to high):");
        _out.WriteLine($"  c0 (constant)  = {coefs[0]:G10}");
        _out.WriteLine($"  c1 (y)         = {coefs[1]:G10}");
        _out.WriteLine($"  c2 (y²)        = {coefs[2]:G10}");
        _out.WriteLine($"  c3 (y³)        = {coefs[3]:G10}");
        _out.WriteLine($"Reduced (× denom guess): try D=18 (path-6 denom):");
        _out.WriteLine($"  18·c0 = {18.0 * coefs[0]:G10} (target: small integer if denom is 18)");
        _out.WriteLine($"  18·c1 = {18.0 * coefs[1]:G10}");
        _out.WriteLine($"  18·c2 = {18.0 * coefs[2]:G10}");
        _out.WriteLine($"  18·c3 = {18.0 * coefs[3]:G10}");

        // Also try common cyclotomic Φ_9-related denominators
        foreach (int d in new[] { 9, 18, 27, 36, 81, 162 })
        {
            _out.WriteLine($"  {d}·c0 = {d * coefs[0]:G10}, {d}·c1 = {d * coefs[1]:G10}, {d}·c2 = {d * coefs[2]:G10}, {d}·c3 = {d * coefs[3]:G10}");
        }
    }

    [Theory]
    [InlineData(5, 2)]    // path-4
    [InlineData(6, 2)]    // path-5
    [InlineData(6, 4)]
    [InlineData(7, 2)]    // path-6
    [InlineData(7, 4)]
    [InlineData(8, 2)]    // path-7
    [InlineData(8, 6)]
    public void ComputeAnatomy_MklDirect_AgreesWithMathNetOracle(int N, int n)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);

        // MklDirect path (the production ComputeAnatomy) via ClassFixture cache.
        var anatomy = _cache.Get(N);
        double? sigmaMklDirect = anatomy.SigmaForBlochIndex(n);
        Assert.NotNull(sigmaMklDirect);

        // MathNet oracle path via reflection (ComputeAnatomyMathNet is private static).
        var mathNetMethod = typeof(C2FullBlockSigmaAnatomy).GetMethod(
            "ComputeAnatomyMathNet",
            System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static);
        Assert.NotNull(mathNetMethod);
        var oracle = (IReadOnlyList<SigmaModeWitness>)mathNetMethod!.Invoke(
            null, new object[] { block, C2FullBlockSigmaAnatomy.DefaultQ })!;
        double? sigmaOracle = null;
        foreach (var w in oracle)
            if (w.BlochIndexN == n) { sigmaOracle = w.Sigma; break; }
        Assert.NotNull(sigmaOracle);

        Assert.True(Math.Abs(sigmaMklDirect!.Value - sigmaOracle!.Value) < 1e-9,
            $"path-{N-1} n={n}: MklDirect={sigmaMklDirect.Value:G10}, MathNet oracle={sigmaOracle.Value:G10}");
    }

    /// <summary>F89 D_k cross-check at well-conditioned stretch (k=25..30): anatomy
    /// inverse-iteration extraction agrees with <see cref="F89UnifiedFaClosedFormClaim.PredictDenominator"/>
    /// within Vandermonde-conditioning floor.
    ///
    /// <para><b>Red signal closed 2026-05-15:</b> the previous k=31, 32 rows that
    /// deviated ~1.5-2·10⁻⁴ against the 1e-4 tolerance are now removed. The red
    /// signal was "something still missing on the route" between extraction and the
    /// empirical D_k formula; the missing piece is the algebraic derivation via the
    /// Chebyshev expansion + orbit-polynomial-reduction pipeline
    /// (<see cref="F89PathPolynomialPipeline"/>), Tier-1-Derived closure. The
    /// pipeline gives bit-exact (P_k, D_k) for arbitrary k; the precision wall at
    /// k=31, 32 is now understood as the extraction *instrument's* conditioning
    /// floor, not a closed-form deficit. See
    /// <c>PredictDenominator_AtAnyK_MatchesPipelineBitExact</c> below for the
    /// authoritative bit-exact verification.</para>
    ///
    /// <para>This test retains the k=25..30 range as a regression check that
    /// <see cref="C2FullBlockSigmaAnatomy.BuildFaOnly"/> + Vandermonde extraction
    /// continues to track the algebraic D_k where conditioning permits. Structural
    /// context: <c>docs/proofs/PROOF_F89_PATH_D_CLOSED_FORM.md</c> "Tier-1-Derived
    /// closure via Chebyshev pipeline".</para>
    /// </summary>
    [Theory]
    [InlineData(25)]   // v₂=0, predicted D = 640000   = 625·2¹⁰
    [InlineData(26)]   // v₂=1, predicted D = 346112   = 169·2¹¹
    [InlineData(27)]   // v₂=0, predicted D = 1492992  = 729·2¹¹
    [InlineData(28)]   // v₂=2, predicted D = 401408   = 49·2¹³
    [InlineData(29)]   // v₂=0, predicted D = 3444736  = 841·2¹²
    [InlineData(30)]   // v₂=1, predicted D = 1843200  = 225·2¹³
    public void PredictDenominator_AtWellConditionedK_MatchesExtractedFromAnatomy(int k)
    {
        int N = k + 1;
        var anatomy = C2FullBlockSigmaAnatomy.BuildFaOnly(C2Block(N));
        double[] rawCoefs = anatomy.ExtractRawPolynomialCoefficients();

        int predictedD = F89UnifiedFaClosedFormClaim.PredictDenominator(k);
        for (int i = 0; i < rawCoefs.Length; i++)
        {
            double scaled = rawCoefs[i] * predictedD;
            double rounded = Math.Round(scaled);
            Assert.True(Math.Abs(scaled - rounded) < 1e-4,
                $"path-{k} coef[{i}]: raw={rawCoefs[i]:G10}, raw·D={scaled:G10}, " +
                $"nearest integer={rounded:G10}, residual={Math.Abs(scaled - rounded):G6}, D={predictedD}");

            _out.WriteLine($"  k={k} coef[{i}] · D({predictedD}) = {scaled:G10} ≈ {(long)rounded}");
        }
    }

    /// <summary>F89 D_k bit-exact verification across the full tabulation range and
    /// past the int-typed boundary, using the native Chebyshev pipeline as the
    /// authoritative source. This replaces the ED-based "red signal" verification
    /// at k=31, 32 (closed 2026-05-15): no inverse iteration, no Vandermonde
    /// conditioning, exact BigInteger equality.</summary>
    [Theory]
    [InlineData(25)]
    [InlineData(31)]
    [InlineData(32)]
    [InlineData(40)]
    [InlineData(46)]
    [InlineData(47)]  // first path past int.MaxValue
    [InlineData(100)]
    public void PredictDenominator_AtAnyK_MatchesPipelineBitExact(int k)
    {
        var pipelineD = F89UnifiedFaClosedFormClaim.ComputePathPolynomialBig(k).Denominator;
        var formulaD = F89UnifiedFaClosedFormClaim.PredictDenominatorBig(k);
        Assert.Equal(pipelineD, formulaD);

        _out.WriteLine($"  k={k}: D_k = {pipelineD} (pipeline ≡ formula, exact)");
    }

    [Theory]
    [InlineData(5)]    // path-4 -> floor(5/2) = 2 F_a modes
    [InlineData(6)]    // path-5 -> 3
    [InlineData(7)]    // path-6 -> 3
    [InlineData(8)]    // path-7 -> 4
    public void BuildFaOnly_ProducesFloorNHalf_FaWitnesses(int N)
    {
        var anatomy = C2FullBlockSigmaAnatomy.BuildFaOnly(C2Block(N));

        Assert.True(anatomy.FaModesOnly);
        Assert.Equal(N / 2, anatomy.SigmaSpectrum.Count);
        Assert.All(anatomy.SigmaSpectrum, w => Assert.True(w.BlochIndexN.HasValue,
            $"BuildFaOnly witness at lambda={w.EigenvalueReal}+{w.EigenvalueImag}i must have a Bloch index"));
        Assert.All(anatomy.SigmaSpectrum, w => Assert.True(w.Sigma > 0,
            $"F_a mode n={w.BlochIndexN} must have sigma > 0; got {w.Sigma}"));
    }

    [Theory]
    [InlineData(3, 2)]   // path-3, N=4
    [InlineData(3, 4)]
    [InlineData(4, 2)]   // path-4, N=5
    [InlineData(4, 4)]
    [InlineData(5, 2)]   // path-5, N=6
    [InlineData(5, 4)]
    [InlineData(5, 6)]
    [InlineData(6, 2)]   // path-6, N=7
    [InlineData(6, 4)]
    [InlineData(6, 6)]
    [InlineData(7, 2)]   // path-7, N=8
    [InlineData(7, 4)]
    [InlineData(7, 6)]
    [InlineData(7, 8)]
    [InlineData(8, 2)]   // path-8, N=9
    [InlineData(8, 4)]
    [InlineData(8, 6)]
    [InlineData(8, 8)]
    [InlineData(9, 2)]   // path-9, N=10
    [InlineData(9, 4)]
    [InlineData(9, 6)]
    [InlineData(9, 8)]
    [InlineData(9, 10)]
    public void BuildFaOnly_InverseIteration_AgreesWithZgeev(int k, int n)
    {
        int N = k + 1;   // C2Block(N) is F89 path-(N-1)
        var block = C2Block(N);

        // zgeev path (full spectrum) — the oracle, via the ClassFixture cache.
        double? sigmaZgeev = _cache.Get(N).SigmaForBlochIndex(n);
        Assert.NotNull(sigmaZgeev);

        // Inverse-iteration path (F_a-only).
        double? sigmaInverseIteration = C2FullBlockSigmaAnatomy.BuildFaOnly(block).SigmaForBlochIndex(n);
        Assert.NotNull(sigmaInverseIteration);

        Assert.True(Math.Abs(sigmaInverseIteration!.Value - sigmaZgeev!.Value) < 1e-9,
            $"path-{k} n={n}: inverse-iteration={sigmaInverseIteration.Value:G10}, " +
            $"zgeev oracle={sigmaZgeev.Value:G10}, " +
            $"diff={Math.Abs(sigmaInverseIteration.Value - sigmaZgeev.Value):G6}");
    }
}

/// <summary>Shared cache of C2FullBlockSigmaAnatomy instances keyed by N, reused
/// across all tests in C2FullBlockSigmaAnatomyTests. Eliminates the per-InlineData
/// redundant eigendecomposition: each unique N triggers Build exactly once across
/// the test class run. Thread-safe (ConcurrentDictionary + Lazy).</summary>
public sealed class C2FullBlockSigmaAnatomyCache
{
    private readonly System.Collections.Concurrent.ConcurrentDictionary<int, Lazy<C2FullBlockSigmaAnatomy>> _byN = new();

    public C2FullBlockSigmaAnatomy Get(int N) =>
        _byN.GetOrAdd(N, n => new Lazy<C2FullBlockSigmaAnatomy>(
            () => C2FullBlockSigmaAnatomy.Build(new CoherenceBlock(N: n, n: 1, gammaZero: 0.05)),
            System.Threading.LazyThreadSafetyMode.ExecutionAndPublication)).Value;
}
