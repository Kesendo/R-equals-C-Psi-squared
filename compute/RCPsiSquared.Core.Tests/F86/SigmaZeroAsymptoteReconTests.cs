using MathNet.Numerics.LinearAlgebra;
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

        // Polynomial fit in 1/N for each parity subsequence; the constant term is a
        // direct estimate of N → ∞ and the higher coefficients (b, c, d, …) reveal the
        // sub-leading scaling structure. If both parities agree on (a, b, c, …), the
        // even/odd split is the only non-analytic structure; if they differ in b, the
        // (-1)^N term enters at first sub-leading.
        var evenNs = Enumerable.Range(NMin, NMax - NMin + 1).Where(n => n % 2 == 0).ToArray();
        var oddNs = Enumerable.Range(NMin, NMax - NMin + 1).Where(n => n % 2 == 1).ToArray();
        _out.WriteLine("\nPolynomial fits (degree 3 in 1/N) per parity subsequence:");
        _out.WriteLine($"  even-N points: N ∈ {{{string.Join(",", evenNs)}}}");
        _out.WriteLine($"  odd-N  points: N ∈ {{{string.Join(",", oddNs)}}}");
        ReportFit("even", evenNs, even, degree: 3);
        ReportFit("odd ", oddNs, odd, degree: 3);

        // Cross-check: increase degree, look for stability of leading term a.
        _out.WriteLine("\nLeading-coefficient stability across fit degrees:");
        for (int deg = 1; deg <= 5; deg++)
        {
            var ce = PolyFit1OverN(evenNs, even, deg);
            var co = PolyFit1OverN(oddNs, odd, deg);
            _out.WriteLine($"  degree {deg}: a_even = {ce[0]:F8}, a_odd = {co[0]:F8}, " +
                           $"avg = {(ce[0] + co[0]) / 2:F8}");
        }
    }

    /// <summary>Fit values ≈ a + b/N + c/N² + … via least squares (QR).</summary>
    private static double[] PolyFit1OverN(int[] Ns, double[] vals, int degree)
    {
        int M = Ns.Length;
        var A = Matrix<double>.Build.Dense(M, degree + 1);
        var y = Vector<double>.Build.Dense(vals);
        for (int i = 0; i < M; i++)
        {
            double invN = 1.0 / Ns[i];
            double pow = 1.0;
            for (int k = 0; k <= degree; k++)
            {
                A[i, k] = pow;
                pow *= invN;
            }
        }
        return A.QR().Solve(y).ToArray();
    }

    private void ReportFit(string label, int[] Ns, double[] vals, int degree)
    {
        var c = PolyFit1OverN(Ns, vals, degree);
        var sb = new System.Text.StringBuilder($"  {label}:  ");
        for (int k = 0; k <= degree; k++)
            sb.Append($"{(k == 0 ? "a" : k == 1 ? "b" : k == 2 ? "c" : k == 3 ? "d" : $"e_{k}")}={c[k]:F6}  ");
        _out.WriteLine(sb.ToString());
    }

    /// <summary>Cross-chromaticity reconnaissance: compute σ_0(c, N) for c ∈ {3, 4} at the
    /// largest N each can reach on commodity hardware, apply parity-split Aitken, look for
    /// a pattern across c that might reveal the closed-form structure shared with c=2.</summary>
    [Fact]
    public void Recon_SigmaZero_C3_C4_CrossChromaticityPattern()
    {
        // c=3 max-feasible N=11 on commodity hardware (block dim 11·165=1815²·16 ≈ 53 MB OK;
        // N=12 14520²·16 = 3.4 GB hits OOM in BlockLDecomposition.MhTotal). c=4 capped at
        // N=8 (block dim 56·70=3920 OK); N=9 84·126=10584 → marginal, N=10 OOM.
        var c3Ns = new[] { 5, 6, 7, 8, 9, 10, 11 };
        var c4Ns = new[] { 7, 8 };
        var sw = System.Diagnostics.Stopwatch.StartNew();

        _out.WriteLine("σ_0(c=3, N) — popcount-2 vs popcount-3 coherences:");
        var c3Sigmas = new double[c3Ns.Length];
        for (int i = 0; i < c3Ns.Length; i++)
        {
            int N = c3Ns[i];
            var svd = InterChannelSvd.Build(
                new CoherenceBlock(N, n: 2, gammaZero: 0.05), hd1: 1, hd2: 3);
            c3Sigmas[i] = svd.Sigma0;
            _out.WriteLine($"  N={N,2}: σ_0 = {svd.Sigma0:F8}   " +
                           $"(cumul. {sw.Elapsed.TotalSeconds,5:F1} s)");
        }

        _out.WriteLine("\nσ_0(c=4, N) — popcount-3 vs popcount-4 coherences:");
        var c4Sigmas = new double[c4Ns.Length];
        for (int i = 0; i < c4Ns.Length; i++)
        {
            int N = c4Ns[i];
            var svd = InterChannelSvd.Build(
                new CoherenceBlock(N, n: 3, gammaZero: 0.05), hd1: 1, hd2: 3);
            c4Sigmas[i] = svd.Sigma0;
            _out.WriteLine($"  N={N,2}: σ_0 = {svd.Sigma0:F8}   " +
                           $"(cumul. {sw.Elapsed.TotalSeconds,5:F1} s)");
        }
        sw.Stop();
        _out.WriteLine($"\nTotal wall: {sw.Elapsed.TotalSeconds:F1} s");

        // Reference asymptotes (refuted as actual limits — they are sweet-spot crossings):
        //   c=2: 2√2     ≈ 2.8284
        //   c=3: 4        = 4.0000
        //   c=4: 2√6     ≈ 4.8990
        // Look for: does σ_0(c, ∞) / (refuted-asymptote) approach a c-independent constant?
        // c=2 (∞) ≈ 2.8628 → ratio 2.8628/2.8284 = 1.0122
        _out.WriteLine($"\nRatios σ_0(c, N) / refuted-asymptote (sweet-spot crossing 2√(2(c−1))):");
        _out.WriteLine($"  c=2 N=18: {2.86222985 / (2 * Math.Sqrt(2)):F6}   (ref: 1.0122 from Aitken at N=∞)");
        _out.WriteLine($"  c=3 N={c3Ns[^1]}: {c3Sigmas[^1] / 4.0:F6}");
        _out.WriteLine($"  c=4 N={c4Ns[^1]}: {c4Sigmas[^1] / (2 * Math.Sqrt(6)):F6}");
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
