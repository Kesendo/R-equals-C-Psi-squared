using System.Diagnostics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Cost probe before any n = 7 β-exotic run (the survey rule: never reuse a runner blindly,
/// check its N-tuning first). The N=5 certificate is cheap for a reason that does NOT scale: it also
/// proves the remainder-R1 gcd, whose resultant R runs against the corner block (p_c+1, p_c+1). That
/// corner has dimension C(N, (N+1)/2+1)² = 25 at N = 5 but 441 at N = 7, so deg_q R jumps from 422 to
/// ≈ 53·441 ≈ 23000 interpolation nodes per prime, times thousands of primes for the Mignotte lift.
///
/// <para>The β-exotic exclusion needs NONE of that: it reads only D(q) = disc_Λ(F_res)(q), whose bound
/// is resDeg·(resDeg−1) ≈ 2756 at N = 7. This probe times the pieces a D-only path would pay for, so
/// the decision to build one rests on measurement rather than on "it is the same call".</para></summary>
public class BetaExoticN7CostProbeTests
{
    private readonly ITestOutputHelper _out;
    public BetaExoticN7CostProbeTests(ITestOutputHelper output) => _out = output;

    /// <summary>Step 1, the cheap half: does the AT factor even reconstruct at n = 7, and how fast?
    /// (At N = 5 the whole per-point path is 0.2-0.7 s; a first n = 7 probe ran past 28 minutes without
    /// finishing, so the pieces are timed separately rather than together.)</summary>
    [Fact]
    [Trait("Category", "BETAN7PROBE")]
    public void N7_AtFactor_IsTimed()
    {
        foreach (bool rOdd in new[] { false, true })
        {
            var sw = Stopwatch.StartNew();
            var at = FoldResultantCertificate.AtFactorAt(7, rOdd, q0: 2);
            _out.WriteLine($"N=7 rOdd={rOdd}: AT degree {at.Length - 1} in {sw.ElapsedMilliseconds} ms");
            Assert.True(at.Length - 1 > 0);
        }
    }

    /// <summary>Step 2: the D-ONLY per-point path, which is what a β-exotic run at n = 7 actually needs.
    /// It is `ExactSample` minus the corner: block charpoly, divide by AT, take the Λ-discriminant. The
    /// corner pencil is 441×441 at n = 7 (25×25 at N = 5) and its ℤ[i] charpoly plus the degree-53×441
    /// resultant are what made the first n = 7 probe run past 28 minutes. Neither is needed here.</summary>
    [Fact]
    [Trait("Category", "BETAN7PROBE")]
    public void N7_DiscOnlyPerPointPath_IsTimed_WithoutTheCorner()
    {
        foreach (int n in new[] { 5, 7 })
        {
            foreach (bool rOdd in new[] { true, false })
            {
                var sw = Stopwatch.StartNew();
                var blk = rOdd
                    ? F89PathKSeDeBlock.BuildTwoTimesROddBlock(2, n)
                    : F89PathKSeDeBlock.BuildTwoTimesSymBlock(2, n);
                long buildMs = sw.ElapsedMilliseconds;

                sw.Restart();
                var chFull = GaussianMatrixCharpoly.Characteristic(blk);
                long charMs = sw.ElapsedMilliseconds;

                sw.Restart();
                var at = FoldResultantCertificate.AtFactorAt(n, rOdd, q0: 2);
                var (res, rem) = GaussianPolynomial.DivMod(chFull, at);
                long divMs = sw.ElapsedMilliseconds;
                Assert.Empty(rem);

                sw.Restart();
                var disc = GaussianPolynomial.Discriminant(res);
                long discMs = sw.ElapsedMilliseconds;

                _out.WriteLine($"n={n} rOdd={rOdd}: block {blk.GetLength(0)}x{blk.GetLength(0)} in {buildMs} ms, " +
                               $"charpoly in {charMs} ms, AT+div (residual degree {res.Length - 1}) in {divMs} ms, " +
                               $"exact disc in {discMs} ms; disc ~ {disc.Re.ToString().Length} decimal digits");
            }
        }
    }

    [Fact]
    [Trait("Category", "BETAN7PROBE")]
    public void N5_SamePieces_ForScaleComparison()
    {
        foreach (bool rOdd in new[] { false, true })
        {
            var sw = Stopwatch.StartNew();
            var (residual, _, _, disc) = FoldResultantCertificate.ExactSample(5, rOdd, q0: 2);
            _out.WriteLine($"N=5 rOdd={rOdd}: residual degree {residual.Length - 1} in {sw.ElapsedMilliseconds} ms; " +
                           $"disc ~ {disc.Re.ToString().Length} decimal digits");
        }
    }
}
