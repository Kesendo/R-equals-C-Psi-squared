using System.Diagnostics;
using RCPsiSquared.Core.F89PathK;
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

    [Fact]
    [Trait("Category", "BETAN7PROBE")]
    public void N7_ExactPerPointPath_IsTimed_BeforeAnyBivariateBuild()
    {
        foreach (bool rOdd in new[] { false, true })
        {
            var sw = Stopwatch.StartNew();
            var at = FoldResultantCertificate.AtFactorAt(7, rOdd, q0: 2);
            long atMs = sw.ElapsedMilliseconds;

            sw.Restart();
            var (residual, _, _, disc) = FoldResultantCertificate.ExactSample(7, rOdd, q0: 2);
            long exMs = sw.ElapsedMilliseconds;

            _out.WriteLine($"N=7 rOdd={rOdd}: AT degree {at.Length - 1} in {atMs} ms; " +
                           $"exact per-point residual degree {residual.Length - 1} in {exMs} ms; " +
                           $"disc value magnitude ~ {disc.Re.ToString().Length} decimal digits (real part)");

            Assert.True(residual.Length - 1 > 0, "the residual must be non-trivial at N = 7");
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
