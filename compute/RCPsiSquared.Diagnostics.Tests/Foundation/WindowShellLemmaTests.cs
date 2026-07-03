using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The WINDOW-COMBINATORICS SHELL LEMMA (PROOF_CODIM1_BY_ADDITIVITY §6, step 2 of the large-N
/// exclusion program): at REAL q the Bendixson rate window confines every block's Re-spectrum by pure
/// popcount combinatorics, N-generally. For the (p,w) block, n_diff = popcount(a⊕b) ranges over
/// [|p−w|, min(p+w, 2N−p−w)] in steps of 2 (parity |p−w| mod 2), so Re spec(p,w)(q ∈ ℝ) ⊆
/// [−2·min(p+w, 2N−p−w), −2·|p−w|]. With the window-edge lemma's N-uniform Re λ_A ∈ (−6,−2):
///
/// <code>
///   λ_A can live only where |p−w| ≤ 2      (the BAND SHELL; margin ≥ Re λ_A + 2|p−w| ≥ Re λ_A + 6 outside)
///   μ = −λ_A−2N only where |p+w−N| ≤ 2     (the HALF-FILLING SHELL, the band shell's s-fold image)
/// </code>
///
/// so at real loci EVERY block outside the two width-5 shells is excluded from {λ_A, μ} for ALL N with the
/// SAME margin Re λ_A + 6 &gt; 0 the corner certificate carries: O(N) survivors instead of O(N²), no
/// certificate, no eigensolve needed. Combined with the fold-lattice quotient (§7, BlockLatticeFoldGroupTests)
/// the survivors live on the fundamental domain's shell strip. The negative control pins where the lemma is
/// SILENT (a shell block whose window contains Re λ_A), which is exactly where the N=5 fold-resultant
/// certificates took over (RemainderR4InteriorExclusionTests).
///
/// <para>Run: <c>dotnet test "compute/RCPsiSquared.Diagnostics.Tests" --filter "Category=WINDOWSHELL"
/// --logger "console;verbosity=detailed"</c>.</para></summary>
public sealed class WindowShellLemmaTests
{
    private readonly ITestOutputHelper _out;
    public WindowShellLemmaTests(ITestOutputHelper o) => _out = o;

    // The two derived real seed loci with their (1,2) defective eigenvalues (proof §6: N=5 q*=0.620878,
    // λ_A=−4.6189; N=7 q*=1.5148, λ_A=−4.885, the pkmono --exact reading).
    private static readonly (int N, double QStar, double LambdaAGuess)[] RealSeeds =
    {
        (5, 0.620878, -4.6189),
        (7, 1.5148, -4.885),
    };

    [Fact]
    [Trait("Category", "WINDOWSHELL")]
    public void NDiffRange_IsPopcountCombinatorics_AllBlocks_N5toN8()
    {
        // n_min = |p−w|, n_max = min(p+w, 2N−p−w), steps of 2: enumerated exactly, no linear algebra.
        for (int n = 5; n <= 8; n++)
            for (int p = 0; p <= n; p++)
                for (int w = 0; w <= n; w++)
                {
                    var kets = Enumerable.Range(0, 1 << n).Where(m => System.Numerics.BitOperations.PopCount((uint)m) == p).ToArray();
                    var bras = Enumerable.Range(0, 1 << n).Where(m => System.Numerics.BitOperations.PopCount((uint)m) == w).ToArray();
                    var seen = new HashSet<int>();
                    foreach (var a in kets) foreach (var b in bras)
                        seen.Add(System.Numerics.BitOperations.PopCount((uint)(a ^ b)));
                    int nMin = Math.Abs(p - w), nMax = Math.Min(p + w, 2 * n - p - w);
                    Assert.Equal(nMin, seen.Min());
                    Assert.Equal(nMax, seen.Max());
                    Assert.All(seen, v => Assert.Equal(nMin % 2, v % 2));
                }
        _out.WriteLine("n_diff range = [|p−w|, min(p+w, 2N−p−w)], parity-stepped: all blocks, N=5..8, exact");
    }

    [Fact]
    [Trait("Category", "WINDOWSHELL")]
    public void RateWindow_HoldsOnEveryBlock_AtRealLoci_N5()
    {
        foreach (var q in new[] { 0.620878, 1.077615 })
        {
            double worstTop = 0, worstBottom = 0;
            for (int p = 0; p <= 5; p++)
                for (int w = 0; w <= 5; w++)
                {
                    var re = Spec(5, p, w, new Complex(q, 0)).Select(z => z.Real).ToArray();
                    int nMin = Math.Abs(p - w), nMax = Math.Min(p + w, 10 - p - w);
                    worstTop = Math.Max(worstTop, re.Max() - (-2.0 * nMin));
                    worstBottom = Math.Max(worstBottom, (-2.0 * nMax) - re.Min());
                }
            _out.WriteLine($"q={q}: window violations top {worstTop:E2}, bottom {worstBottom:E2} (all 36 blocks)");
            Assert.True(worstTop < 1e-9 && worstBottom < 1e-9,
                $"Bendixson must confine every block at real q={q}; top {worstTop:E2}, bottom {worstBottom:E2}");
        }
    }

    [Fact]
    [Trait("Category", "WINDOWSHELL")]
    public void OutsideTheTwoShells_EveryBlockExcludes_LambdaA_AndMu_N5andN7()
    {
        foreach (var (n, qStar, guess) in RealSeeds)
        {
            var q = new Complex(qStar, 0);
            var lamA = Spec(n, 1, 2, q).OrderBy(z => Math.Abs(z.Real - guess) + Math.Abs(z.Imaginary)).First();
            double mu = -lamA.Real - 2 * n;
            double uniformMargin = lamA.Real + 6;
            Assert.True(uniformMargin > 0, "the window-edge lemma's Re λ_A > −6 must hold at the seed");

            int checkedBlocks = 0;
            double worstSlackA = double.PositiveInfinity, worstSlackMu = double.PositiveInfinity;
            for (int p = 0; p <= n; p++)
                for (int w = 0; w <= n; w++)
                {
                    bool outsideBand = Math.Abs(p - w) >= 3;
                    bool outsideHalf = Math.Abs(p + w - n) >= 3;
                    if (!outsideBand && !outsideHalf) continue;
                    var spec = Spec(n, p, w, q);
                    checkedBlocks++;
                    if (outsideBand)
                    {
                        double margin = lamA.Real + 2.0 * Math.Abs(p - w);        // ≥ uniformMargin ≥ Re λ_A + 6
                        double d = spec.Min(z => (z - lamA).Magnitude);
                        worstSlackA = Math.Min(worstSlackA, d - margin);
                        Assert.True(d >= margin - 1e-9,
                            $"N={n} ({p},{w}) |p−w|≥3 must exclude λ_A with the combinatorial margin; dist {d:F4} < {margin:F4}");
                    }
                    if (outsideHalf)
                    {
                        double margin = uniformMargin;                            // −2(N−3) − Re μ = Re λ_A + 6
                        double d = spec.Min(z => Math.Abs(z.Real - mu));
                        worstSlackMu = Math.Min(worstSlackMu, d - margin);
                        Assert.True(d >= margin - 1e-9,
                            $"N={n} ({p},{w}) |p+w−N|≥3 must exclude μ with margin Re λ_A + 6; dist {d:F4} < {margin:F4}");
                    }
                }
            _out.WriteLine($"N={n} q*={qStar}: λ_A={lamA.Real:F4}, μ={mu:F4}, uniform margin {uniformMargin:F4}; " +
                           $"{checkedBlocks} outside-shell blocks excluded (worst slack λ_A {worstSlackA:F3}, μ {worstSlackMu:F3})");
        }
    }

    [Fact]
    [Trait("Category", "WINDOWSHELL")]
    public void NegativeControl_InsideTheShell_TheWindowIsSilent()
    {
        // (1,1) at the second N=5 real locus: λ_A = −3.7917 sits INSIDE the (1,1) window [−4, 0], so the
        // combinatorial lemma cannot exclude it; that residual shell is exactly what the fold-resultant
        // certificate pair closed (RemainderR4InteriorExclusionTests). The lemma's silence is by design.
        var q = new Complex(1.077615, 0);
        var lamA = Spec(5, 1, 2, q).OrderBy(z => Math.Abs(z.Real + 3.7917) + Math.Abs(z.Imaginary)).First();
        Assert.InRange(lamA.Real, -2.0 * 2, -2.0 * 0);          // inside the (1,1) window [−2·n_max, −2·n_min] = [−4, 0]
        _out.WriteLine($"locus 2: Re λ_A = {lamA.Real:F4} ∈ [−4, 0], the (1,1) window: the shell needs the certificate");
    }

    private static Complex[] Spec(int n, int wk, int wb, Complex q) =>
        Matrix<Complex>.Build.DenseOfArray(WeightCoherenceBlock.Build(n, wk, wb, q)).Evd().EigenValues.ToArray();
}
