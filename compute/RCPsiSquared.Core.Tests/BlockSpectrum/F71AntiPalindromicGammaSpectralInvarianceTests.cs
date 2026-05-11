using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

/// <summary>Tests for <see cref="F71AntiPalindromicGammaSpectralInvariance"/>: the chain
/// XY+Z-dephasing Liouvillian F71-REFINED DIAGONAL-BLOCK eigenvalue multiset is invariant
/// under any γ-distribution satisfying γ_l + γ_{N-1-l} = 2·γ_avg (anti-palindromic around
/// its mean).
///
/// <para><b>Scope.</b> The invariance applies to the diagonal sub-blocks of <c>Q^T L Q</c>
/// only (computed by <see cref="F71MirrorBlockRefinement.ComputeSpectrumPerBlock"/>). The
/// full L spectrum (computed by <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/>)
/// generally DIFFERS across anti-palindromic γ-profiles because the (F71-even ↔ F71-odd)
/// cross blocks are nonzero when γ breaks F71. We test BOTH facts: diagonal-block
/// invariance + full-L deviation.</para>
///
/// <para>Strictly weaker than F71 symmetry only at uniform γ (the unique overlap of
/// F71-palindromy and anti-palindromy); strictly stronger than F1 (Σγ_l invariant). For
/// odd N the middle site is forced to equal γ_avg by self-pairing.</para></summary>
public class F71AntiPalindromicGammaSpectralInvarianceTests
{
    // ----------------------------------------------------------------------
    // Claim metadata
    // ----------------------------------------------------------------------

    [Fact]
    public void Tier_IsTier1Candidate()
    {
        var claim = new F71AntiPalindromicGammaSpectralInvariance(
            new JointPopcountSectors(),
            new F71MirrorBlockRefinement(new JointPopcountSectors()));
        Assert.Equal(Tier.Tier1Candidate, claim.Tier);
    }

    [Fact]
    public void DisplayName_IsNonEmpty()
    {
        var claim = new F71AntiPalindromicGammaSpectralInvariance(
            new JointPopcountSectors(),
            new F71MirrorBlockRefinement(new JointPopcountSectors()));
        Assert.False(string.IsNullOrWhiteSpace(claim.DisplayName));
    }

    [Fact]
    public void Summary_IsNonEmpty()
    {
        var claim = new F71AntiPalindromicGammaSpectralInvariance(
            new JointPopcountSectors(),
            new F71MirrorBlockRefinement(new JointPopcountSectors()));
        Assert.False(string.IsNullOrWhiteSpace(claim.Summary));
    }

    [Fact]
    public void Constructor_NullSectors_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F71AntiPalindromicGammaSpectralInvariance(
                null!, new F71MirrorBlockRefinement(new JointPopcountSectors())));
    }

    [Fact]
    public void Constructor_NullF71_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F71AntiPalindromicGammaSpectralInvariance(new JointPopcountSectors(), null!));
    }

    // ----------------------------------------------------------------------
    // AntiPalindromicDeviation: closed-form correctness on known γ-distributions
    // ----------------------------------------------------------------------

    [Fact]
    public void AntiPalindromicDeviation_NullList_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            F71AntiPalindromicGammaSpectralInvariance.AntiPalindromicDeviation(null!));
    }

    [Fact]
    public void AntiPalindromicDeviation_Uniform_IsZero()
    {
        var gammas = new[] { 0.45, 0.45, 0.45, 0.45, 0.45, 0.45 };
        double dev = F71AntiPalindromicGammaSpectralInvariance.AntiPalindromicDeviation(gammas);
        Assert.True(dev < 1e-12, $"Uniform γ: expected anti-palindromic deviation ≈ 0, got {dev:E3}.");
    }

    [Fact]
    public void AntiPalindromicDeviation_PalindromicNonFlat_IsPositive()
    {
        // [0.3, 0.4, 0.5, 0.4, 0.3] — F71-palindromic (γ_l = γ_{N-1-l}) but γ_avg = 0.38.
        // pairs (0.3+0.3)=0.6, (0.4+0.4)=0.8, (0.5+0.5)=1.0, (0.4+0.4)=0.8, (0.3+0.3)=0.6.
        // 2·γ_avg = 0.76; pair-sum-vs-2avg deviations²:
        //   2·(0.6−0.76)² + 2·(0.8−0.76)² + (1.0−0.76)² = 2·0.0256 + 2·0.0016 + 0.0576 = 0.1120.
        // Demonstrates: F71-palindromy is NOT a subset of anti-palindromy unless γ is flat.
        var palindromicNonFlat = new[] { 0.3, 0.4, 0.5, 0.4, 0.3 };
        double dev = F71AntiPalindromicGammaSpectralInvariance.AntiPalindromicDeviation(palindromicNonFlat);
        Assert.True(dev > 0.3,
            $"Palindromic-but-not-flat γ: expected deviation > 0.3 (palindromy ⊄ anti-palindromy), got {dev:E3}.");
    }

    [Fact]
    public void AntiPalindromicDeviation_AntiPalindromicMonotonic_N6_IsZero()
    {
        var gammas = new[] { 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 };
        double dev = F71AntiPalindromicGammaSpectralInvariance.AntiPalindromicDeviation(gammas);
        Assert.True(dev < 1e-12, $"Monotonic anti-palindromic γ: expected deviation ≈ 0, got {dev:E3}.");
    }

    [Fact]
    public void AntiPalindromicDeviation_AntiPalindromicNonMonotonic_N6_IsZero()
    {
        var gammas = new[] { 0.3, 0.5, 0.4, 0.5, 0.4, 0.6 };
        double dev = F71AntiPalindromicGammaSpectralInvariance.AntiPalindromicDeviation(gammas);
        Assert.True(dev < 1e-12, $"Non-monotonic anti-palindromic γ: expected deviation ≈ 0, got {dev:E3}.");
    }

    [Fact]
    public void AntiPalindromicDeviation_PermutedNonAntiPalindromic_N6_IsPositive()
    {
        // [0.7, 0.2, 0.5, 0.3, 0.6, 0.4] — same γ-multiset as monotonic, γ_avg = 0.45,
        // pairs (0.7+0.4)=1.1, (0.2+0.6)=0.8, (0.5+0.3)=0.8, (0.3+0.5)=0.8, (0.6+0.2)=0.8, (0.4+0.7)=1.1.
        // 2·γ_avg = 0.9; deviations² = 2·(0.2)² + 4·(−0.1)² = 2·0.04 + 4·0.01 = 0.12.
        // sqrt(0.12) ≈ 0.3464.
        var gammas = new[] { 0.7, 0.2, 0.5, 0.3, 0.6, 0.4 };
        double dev = F71AntiPalindromicGammaSpectralInvariance.AntiPalindromicDeviation(gammas);
        Assert.True(dev > 0.3, $"Permuted non-anti-palindromic γ: expected deviation > 0.3, got {dev:E3}.");
    }

    [Fact]
    public void AntiPalindromicDeviation_Concentrated_N6_IsLarge()
    {
        // [0.1, 0.1, 0.1, 0.1, 0.1, 2.2] — γ_avg = 0.45, 2·γ_avg = 0.9.
        // pairs (0.1+2.2)=2.3 ×2 + (0.1+0.1)=0.2 ×4. deviations² = 2·(1.4)² + 4·(−0.7)² = 2·1.96 + 4·0.49 = 5.88.
        // sqrt(5.88) ≈ 2.4248.
        var gammas = new[] { 0.1, 0.1, 0.1, 0.1, 0.1, 2.2 };
        double dev = F71AntiPalindromicGammaSpectralInvariance.AntiPalindromicDeviation(gammas);
        Assert.True(dev > 2.0, $"Concentrated γ: expected deviation > 2.0, got {dev:E3}.");
    }

    // ----------------------------------------------------------------------
    // IsAntiPalindromic: boolean predicate sanity checks
    // ----------------------------------------------------------------------

    [Fact]
    public void IsAntiPalindromic_Uniform_True()
    {
        var gammas = new[] { 0.5, 0.5, 0.5, 0.5 };
        Assert.True(F71AntiPalindromicGammaSpectralInvariance.IsAntiPalindromic(gammas));
    }

    [Fact]
    public void IsAntiPalindromic_MonotonicAntiPalindromic_N6_True()
    {
        var gammas = new[] { 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 };
        Assert.True(F71AntiPalindromicGammaSpectralInvariance.IsAntiPalindromic(gammas));
    }

    [Fact]
    public void IsAntiPalindromic_NonMonotonicAntiPalindromic_N6_True()
    {
        var gammas = new[] { 0.3, 0.5, 0.4, 0.5, 0.4, 0.6 };
        Assert.True(F71AntiPalindromicGammaSpectralInvariance.IsAntiPalindromic(gammas));
    }

    [Fact]
    public void IsAntiPalindromic_Permuted_N6_False()
    {
        var gammas = new[] { 0.7, 0.2, 0.5, 0.3, 0.6, 0.4 };
        Assert.False(F71AntiPalindromicGammaSpectralInvariance.IsAntiPalindromic(gammas));
    }

    [Fact]
    public void IsAntiPalindromic_OddN_MiddleSiteMustEqualMean()
    {
        // [0.3, 0.4, 0.5, 0.6, 0.7] — γ_avg = 0.5, middle site = 0.5 ✓, all pair sums = 1.0 = 2·γ_avg.
        var ok = new[] { 0.3, 0.4, 0.5, 0.6, 0.7 };
        Assert.True(F71AntiPalindromicGammaSpectralInvariance.IsAntiPalindromic(ok));

        // [0.3, 0.4, 0.6, 0.6, 0.7] — γ_avg = 2.6/5 = 0.52, middle = 0.6 ≠ 0.52.
        // pairs (0.3+0.7)=1.0, (0.4+0.6)=1.0, (0.6+0.6)=1.2, (0.6+0.4)=1.0, (0.7+0.3)=1.0;
        // 2·γ_avg = 1.04; deviations² = 4·(−0.04)² + (0.16)² = 0.0064 + 0.0256 = 0.032; sqrt ≈ 0.179.
        var bad = new[] { 0.3, 0.4, 0.6, 0.6, 0.7 };
        Assert.False(F71AntiPalindromicGammaSpectralInvariance.IsAntiPalindromic(bad));
    }

    // ----------------------------------------------------------------------
    // F71-REFINED DIAGONAL-BLOCK spectral invariance: anti-palindromic γ ⇒ same multiset
    // ----------------------------------------------------------------------

    [Fact]
    public void F71RefinedSpectrum_AntiPalindromicMonotonic_N4_MatchesUniformAvg()
    {
        // [0.3, 0.4, 0.5, 0.6] — γ_avg = 0.45, pairs all 0.9.
        const int N = 4;
        var antiPal = new[] { 0.3, 0.4, 0.5, 0.6 };
        var uniform = Enumerable.Repeat(0.45, N).ToArray();
        AssertF71RefinedSpectraEqual(N, uniform, antiPal, binWidth: 1e-7);
    }

    [Fact]
    public void F71RefinedSpectrum_AntiPalindromicMonotonic_N5_MatchesUniformAvg()
    {
        // [0.25, 0.35, 0.45, 0.55, 0.65] — γ_avg = 0.45, middle = 0.45 ✓, pairs all 0.9.
        const int N = 5;
        var antiPal = new[] { 0.25, 0.35, 0.45, 0.55, 0.65 };
        var uniform = Enumerable.Repeat(0.45, N).ToArray();
        AssertF71RefinedSpectraEqual(N, uniform, antiPal, binWidth: 1e-7);
    }

    [Fact]
    public void F71RefinedSpectrum_AntiPalindromicMonotonic_N6_MatchesUniformAvg()
    {
        // [0.2, 0.3, 0.4, 0.5, 0.6, 0.7] — γ_avg = 0.45, pairs all 0.9.
        const int N = 6;
        var antiPal = new[] { 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 };
        var uniform = Enumerable.Repeat(0.45, N).ToArray();
        AssertF71RefinedSpectraEqual(N, uniform, antiPal, binWidth: 1e-7);
    }

    [Fact]
    public void F71RefinedSpectrum_AntiPalindromicNonMonotonic_N6_MatchesUniformAvg()
    {
        // [0.3, 0.5, 0.4, 0.5, 0.4, 0.6] — γ_avg = 0.45, pairs all 0.9.
        const int N = 6;
        var antiPal = new[] { 0.3, 0.5, 0.4, 0.5, 0.4, 0.6 };
        var uniform = Enumerable.Repeat(0.45, N).ToArray();
        AssertF71RefinedSpectraEqual(N, uniform, antiPal, binWidth: 1e-7);
    }

    // ----------------------------------------------------------------------
    // F71-REFINED DIAGONAL-BLOCK spectral breaking: non-anti-palindromic γ ⇒ differs
    // ----------------------------------------------------------------------

    [Fact]
    public void F71RefinedSpectrum_PermutedNonAntiPalindromic_N6_DiffersFromUniformAvg()
    {
        // [0.7, 0.2, 0.5, 0.3, 0.6, 0.4] — same γ-multiset as monotonic but pairs not constant.
        const int N = 6;
        var permuted = new[] { 0.7, 0.2, 0.5, 0.3, 0.6, 0.4 };
        var uniform = Enumerable.Repeat(0.45, N).ToArray();
        int diffBins = CountF71RefinedHistogramMismatchBins(N, uniform, permuted, binWidth: 1e-3);
        Assert.True(diffBins > 0,
            $"N={N}: permuted non-anti-palindromic γ produced 0 mismatched bins; expected > 0 (significant breaking).");
    }

    [Fact]
    public void F71RefinedSpectrum_Concentrated_N6_DiffersFromUniformAvg()
    {
        // [0.1, 0.1, 0.1, 0.1, 0.1, 2.2] — Σγ = 2.7 but heavily skewed.
        const int N = 6;
        var concentrated = new[] { 0.1, 0.1, 0.1, 0.1, 0.1, 2.2 };
        var uniform = Enumerable.Repeat(0.45, N).ToArray();
        int diffBins = CountF71RefinedHistogramMismatchBins(N, uniform, concentrated, binWidth: 1e-3);
        Assert.True(diffBins > 50,
            $"N={N}: concentrated γ produced {diffBins} mismatched bins; expected > 50 (large deformation).");
    }

    // ----------------------------------------------------------------------
    // Full-L spectrum DIFFERS even under anti-palindromic γ — invariance is
    // F71-projected only, not full-L.
    // ----------------------------------------------------------------------

    [Fact]
    public void FullLSpectrum_AntiPalindromicMonotonic_N5_DiffersFromUniformAvg()
    {
        // The diagonal-block invariance is NOT a full-L invariance: the F71 cross blocks
        // (carrying the F71-asymmetry, see InhomogeneousGammaF71BreakingWitness) perturb
        // the full L spectrum away from the diagonal-only multiset.
        const int N = 5;
        var antiPal = new[] { 0.25, 0.35, 0.45, 0.55, 0.65 };
        var uniform = Enumerable.Repeat(0.45, N).ToArray();
        int diffBins = CountFullLHistogramMismatchBins(N, uniform, antiPal, binWidth: 1e-3);
        Assert.True(diffBins > 0,
            $"N={N}: anti-palindromic γ but full-L spectra are equal? Got {diffBins} mismatched bins; expected > 0 (cross-block perturbation should lift degeneracies).");
    }

    [Fact]
    public void FullLSpectrum_AntiPalindromicMonotonic_N6_DiffersFromUniformAvg()
    {
        const int N = 6;
        var antiPal = new[] { 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 };
        var uniform = Enumerable.Repeat(0.45, N).ToArray();
        int diffBins = CountFullLHistogramMismatchBins(N, uniform, antiPal, binWidth: 1e-3);
        Assert.True(diffBins > 0,
            $"N={N}: anti-palindromic γ but full-L spectra are equal? Got {diffBins} mismatched bins; expected > 0 (cross-block perturbation should lift degeneracies).");
    }

    // ----------------------------------------------------------------------
    // F71 BROKEN in operator (off-block Frobenius ≠ 0) but F71-refined diagonal spectrum
    // STAYS invariant — the breaking lives entirely in the cross blocks.
    // ----------------------------------------------------------------------

    [Fact]
    public void AntiPalindromicMonotonic_N6_BreaksF71InOperator_ButPreservesF71RefinedDiagSpectrum()
    {
        const int N = 6;
        var antiPal = new[] { 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 };

        // F71-asymmetry norm: monotonic profile has γ_l ≠ γ_{N-1-l} → asymmetric.
        double f71Asym = InhomogeneousGammaF71BreakingWitness.F71AsymmetryNorm(antiPal);
        Assert.True(f71Asym > 0.1,
            $"Monotonic anti-palindromic γ must be F71-asymmetric in operator form; got asymmetry norm = {f71Asym:E3}.");

        // Anti-palindromic deviation = 0.
        double antiPalDev = F71AntiPalindromicGammaSpectralInvariance.AntiPalindromicDeviation(antiPal);
        Assert.True(antiPalDev < 1e-12,
            $"Monotonic anti-palindromic γ must satisfy anti-palindromy; got deviation = {antiPalDev:E3}.");

        // F71-refined diagonal-block spectrum is invariant despite operator-level F71 breaking.
        var uniform = Enumerable.Repeat(0.45, N).ToArray();
        AssertF71RefinedSpectraEqual(N, uniform, antiPal, binWidth: 1e-7);
    }

    // ----------------------------------------------------------------------
    // Helpers
    // ----------------------------------------------------------------------

    /// <summary>Asserts that the F71-refined diagonal-block eigenvalue multisets of L
    /// (computed by <see cref="F71MirrorBlockRefinement.ComputeSpectrumPerBlock"/>) for
    /// two γ-distributions are equal as histograms binned at <paramref name="binWidth"/>.</summary>
    private static void AssertF71RefinedSpectraEqual(int N, double[] gammaA, double[] gammaB, double binWidth)
    {
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var spectrumA = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, gammaA, N);
        var spectrumB = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, gammaB, N);

        Assert.Equal(spectrumA.Length, spectrumB.Length);
        if (!HistogramsEqual(spectrumA, spectrumB, binWidth, out string? diff))
            Assert.Fail($"N={N}: F71-refined diagonal-block spectra differ as histograms (binWidth={binWidth:E1}). {diff}");
    }

    /// <summary>Counts the number of histogram bins where two F71-refined diagonal-block
    /// spectra differ in count. Returns 0 iff the multisets are equal at this bin width.</summary>
    private static int CountF71RefinedHistogramMismatchBins(int N, double[] gammaA, double[] gammaB, double binWidth)
    {
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var spectrumA = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, gammaA, N);
        var spectrumB = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, gammaB, N);

        Assert.Equal(spectrumA.Length, spectrumB.Length);
        return CountHistogramMismatchBins(spectrumA, spectrumB, binWidth);
    }

    /// <summary>Counts mismatched histogram bins for the FULL L spectrum (joint-popcount
    /// block-eig path, no F71 truncation) between two γ-distributions.</summary>
    private static int CountFullLHistogramMismatchBins(int N, double[] gammaA, double[] gammaB, double binWidth)
    {
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var spectrumA = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaA, N);
        var spectrumB = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaB, N);

        Assert.Equal(spectrumA.Length, spectrumB.Length);
        return CountHistogramMismatchBins(spectrumA, spectrumB, binWidth);
    }

    /// <summary>Bin-based multiset equality: snap each complex eigenvalue to a (binRe, binIm)
    /// lattice with bin width = <paramref name="binWidth"/>, then compare counts. Returns
    /// true iff every bin has identical counts in both. Avoids both sort-key Im-flip
    /// pathologies and greedy-matching degeneracy chains.</summary>
    private static bool HistogramsEqual(IReadOnlyList<Complex> a, IReadOnlyList<Complex> b, double binWidth, out string? diff)
    {
        diff = null;
        if (a.Count != b.Count) { diff = $"length: {a.Count} vs {b.Count}"; return false; }
        var histA = BuildHistogram(a, binWidth);
        var histB = BuildHistogram(b, binWidth);
        var allKeys = new HashSet<(long, long)>(histA.Keys);
        allKeys.UnionWith(histB.Keys);
        var diffs = new List<string>();
        foreach (var key in allKeys)
        {
            int ca = histA.GetValueOrDefault(key);
            int cb = histB.GetValueOrDefault(key);
            if (ca != cb)
            {
                double rE = key.Item1 * binWidth;
                double iE = key.Item2 * binWidth;
                diffs.Add($"bin (Re≈{rE:G6}, Im≈{iE:G6}): A_count={ca}, B_count={cb}");
            }
        }
        if (diffs.Count == 0) return true;
        diff = string.Join("; ", diffs.Take(5)) + (diffs.Count > 5 ? $"; ... +{diffs.Count - 5} more" : "");
        return false;
    }

    private static int CountHistogramMismatchBins(IReadOnlyList<Complex> a, IReadOnlyList<Complex> b, double binWidth)
    {
        var histA = BuildHistogram(a, binWidth);
        var histB = BuildHistogram(b, binWidth);
        var allKeys = new HashSet<(long, long)>(histA.Keys);
        allKeys.UnionWith(histB.Keys);
        int mismatched = 0;
        foreach (var key in allKeys)
        {
            if (histA.GetValueOrDefault(key) != histB.GetValueOrDefault(key)) mismatched++;
        }
        return mismatched;
    }

    private static Dictionary<(long, long), int> BuildHistogram(IReadOnlyList<Complex> spectrum, double binWidth)
    {
        var hist = new Dictionary<(long, long), int>();
        foreach (var z in spectrum)
        {
            var key = ((long)Math.Round(z.Real / binWidth), (long)Math.Round(z.Imaginary / binWidth));
            hist[key] = hist.GetValueOrDefault(key) + 1;
        }
        return hist;
    }
}
