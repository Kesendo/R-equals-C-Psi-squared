using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>F86 Item 5 reconnaissance: extend the σ_0(c=2, N) data table beyond N=11
/// (the current upper bound in <see cref="RCPsiSquared.Core.F86.SigmaZeroChromaticityScaling"/>)
/// and apply Aitken/Richardson extrapolation to sharpen the N → ∞ asymptote estimate.
/// Baseline from existing data: monotone growth 2.83 → 2.86 over N=7..11, Aitken gives
/// ~2.85..2.89. Goal here: more data points + iterated Aitken to narrow the band.</summary>
public class SigmaZeroAsymptoteReconTests
{
    private readonly ITestOutputHelper _out;
    public SigmaZeroAsymptoteReconTests(ITestOutputHelper output) => _out = output;

    /// <summary>Compute σ_0(c=2, N) for N up to <see cref="NMax"/> and apply iterated
    /// Aitken Δ² extrapolation. Writes a table for inspection; asserts the sequence is
    /// monotonically increasing (a structural property already noted in
    /// <see cref="RCPsiSquared.Core.F86.F86OpenQuestions"/> Item 5).</summary>
    private const int NMax = 18;
    private const int NMin = 5;
    private const double GammaZero = 0.05;

    [Fact]
    public void Recon_SigmaZero_C2_GammaIndependence_Confirms_HOnly_Origin()
    {
        // V_inter = P_HD1† · M_H_total · P_HD2 uses only the Hamiltonian part of the
        // Liouvillian (M_H, not D). σ_0 should be γ-independent — verify by running
        // two γ values at the same N and asserting bit-equality (within FP noise).
        const int N = 10;
        var sIdx1 = InterChannelSvd.Build(new CoherenceBlock(N, n: 1, gammaZero: 0.01), 1, 3).Sigma0;
        var sIdx2 = InterChannelSvd.Build(new CoherenceBlock(N, n: 1, gammaZero: 0.50), 1, 3).Sigma0;
        var sIdx3 = InterChannelSvd.Build(new CoherenceBlock(N, n: 1, gammaZero: 5.00), 1, 3).Sigma0;
        _out.WriteLine($"σ_0(c=2, N=10) at γ=0.01: {sIdx1:F12}");
        _out.WriteLine($"σ_0(c=2, N=10) at γ=0.50: {sIdx2:F12}");
        _out.WriteLine($"σ_0(c=2, N=10) at γ=5.00: {sIdx3:F12}");
        Assert.Equal(sIdx1, sIdx2, precision: 12);
        Assert.Equal(sIdx1, sIdx3, precision: 12);
    }

    [Fact]
    public void Recon_SigmaZero_C2_ExtendedNRange_MonotoneAndAitkenConverges()
    {
        var sigmas = new double[NMax - NMin + 1];
        var sw = System.Diagnostics.Stopwatch.StartNew();
        for (int N = NMin; N <= NMax; N++)
        {
            var block = new CoherenceBlock(N, n: 1, gammaZero: GammaZero);
            var svd = InterChannelSvd.Build(block, hd1: 1, hd2: 3);
            sigmas[N - NMin] = svd.Sigma0;
            _out.WriteLine($"  σ_0(c=2, N={N,2}) = {svd.Sigma0:F10}  " +
                           $"(elapsed {sw.Elapsed.TotalSeconds,5:F1} s, " +
                           $"c=2 block dim {N * N * (N - 1) / 2})");
        }
        sw.Stop();
        _out.WriteLine($"Total wall time: {sw.Elapsed.TotalSeconds:F1} s");

        // Structural property: monotone increasing in N (per Item 5 metadata).
        for (int i = 1; i < sigmas.Length; i++)
            Assert.True(sigmas[i] > sigmas[i - 1],
                $"Monotonicity broken at N={NMin + i}: {sigmas[i - 1]:F10} → {sigmas[i]:F10}");

        // Aitken Δ² extrapolation: for a sequence a_n converging to L with leading error
        // proportional to ρ^n, the transform â_n = a_{n+1} − (a_{n+1}−a_n)²/(a_{n+1}−2a_n+a_{n-1})
        // accelerates convergence to L. Iterate for stronger acceleration.
        _out.WriteLine("\nAitken iterates:");
        var aitken1 = AitkenTransform(sigmas);
        var aitken2 = AitkenTransform(aitken1);
        var aitken3 = AitkenTransform(aitken2);

        PrintIterate("a_n      ", sigmas, NMin);
        PrintIterate("â_n  (1st)", aitken1, NMin + 1);
        PrintIterate("â̂_n (2nd)", aitken2, NMin + 2);
        PrintIterate("â̂̂_n (3rd)", aitken3, NMin + 3);

        // Richardson: assume a_N = a_∞ + b/N + O(1/N²); solve b, a_∞ from any pair (a_N, a_M).
        _out.WriteLine("\nRichardson (1/N model) using last two points:");
        double aN1 = sigmas[^2], aN2 = sigmas[^1];
        int N1 = NMax - 1, N2 = NMax;
        double bFit = (aN2 - aN1) * N1 * N2 / (double)(N1 - N2);
        double linfFit = aN2 - bFit / N2;
        _out.WriteLine($"  a_∞ ≈ {linfFit:F6}, b ≈ {bFit:F6}  " +
                       $"(model: a_N = a_∞ + b/N)");

        // Even/odd-N parity in the convergence: successive Δ's halve in pairs, suggesting
        // even-N and odd-N converge separately. Aitken on the mixed sequence oscillates;
        // splitting cleans the convergence.
        _out.WriteLine("\nEven-N / odd-N split (parity in convergence rate):");
        var even = sigmas.Where((_, i) => (NMin + i) % 2 == 0).ToArray();
        var odd = sigmas.Where((_, i) => (NMin + i) % 2 == 1).ToArray();
        PrintIterate("even-N    ", even, NMin + (NMin % 2 == 0 ? 0 : 1));
        PrintIterate("odd-N     ", odd, NMin + (NMin % 2 == 1 ? 0 : 1));

        var evenAitk = AitkenTransform(even);
        var oddAitk = AitkenTransform(odd);
        PrintIterate("even Aitk ", evenAitk, NMin + (NMin % 2 == 0 ? 0 : 1) + 2);
        PrintIterate("odd  Aitk ", oddAitk, NMin + (NMin % 2 == 1 ? 0 : 1) + 2);

        // Quick estimate of N → ∞: average of last even-Aitken and last odd-Aitken.
        double avgLimit = (evenAitk[^1] + oddAitk[^1]) / 2.0;
        _out.WriteLine($"\nAvg of last (even-Aitken, odd-Aitken) = {avgLimit:F6}  " +
                       $"— current best numerical estimate of σ_0(c=2, N → ∞).");
    }

    private static double[] AitkenTransform(double[] a)
    {
        if (a.Length < 3) return Array.Empty<double>();
        var b = new double[a.Length - 2];
        for (int i = 0; i < b.Length; i++)
        {
            double d2 = a[i + 2] - 2 * a[i + 1] + a[i];
            b[i] = Math.Abs(d2) < 1e-15
                ? a[i + 2]
                : a[i + 2] - (a[i + 2] - a[i + 1]) * (a[i + 2] - a[i + 1]) / d2;
        }
        return b;
    }

    private void PrintIterate(string label, double[] arr, int startN)
    {
        var sb = new System.Text.StringBuilder($"  {label}:  ");
        for (int i = 0; i < arr.Length; i++)
            sb.Append($"N={startN + i,2}: {arr[i]:F8}  ");
        _out.WriteLine(sb.ToString());
    }
}
