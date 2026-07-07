using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>THE F86a BLIND SPOT, SHOWN FROM BELOW. F86a-retraction (2026-06-21) concluded the full (1,2)
/// coherence block has NO real-axis defective exceptional point (eigenvalues stay simple, a finite-Petermann
/// off-axis shadow only). F89 (2026-07-04) PROVES a real defective seed exists on the (1,2) block at every odd
/// N (Kato simple-zero of the discriminant; census-confirmed to N=11; β-exotic scoping reads p≈0.5). Same
/// block, same real q axis (GLOSSARY §q-and-Q: ‖L_F86(J) − L_F89(J/2)‖ = 0). The two collide, and F89's
/// theorem wins. This test shows WHY F86a's scan missed what F89 proves.
///
/// <para><b>The mechanism: a detection window far narrower than the scan grid.</b> A defective √-EP splits its
/// coalescing pair by ~√|q − q*|. Any "sit at a q and characterize" reader (F86a's max-Petermann-K + Riesz
/// trio) can only SEE the coalescence when it sits INSIDE a narrow window |q − q*| ≲ ε where the split is
/// small. F86a scanned Q ∈ [0.5, 4.0] on ~121 points (ΔQ ≈ 0.029), so its nearest grid point sits up to
/// ~0.0146 from any q*. This test measures the coalescing pair's gap at that offset directly: it is ~1e-1
/// there (two plainly-separate eigenvalues), and only collapses toward 0 within ~1e-3 of q*. So the scan
/// reported "no real-axis defective EP" at every grid point, and the sharp Petermann swing at the true q*
/// (the cusp) was dismissed as a grid artifact.</para>
///
/// <para><b>Why F89 sees it: a grid-robust detector.</b> F89 does NOT sit-and-characterize as its primary
/// finder. It counts REAL residual roots (<see cref="PathKMonodromyScout.FindRealDefectiveByCountChange"/>):
/// across q* two real strands merge and leave the axis, so the integer real-root count jumps by 2 — a
/// TOPOLOGICAL signal visible at ANY grid resolution, immune to how far the grid sits from q*. It flags the
/// transition on a coarse grid, THEN bisects onto q* (~1e-7) and only there characterizes. The window problem
/// that blinds a sit-and-characterize scan is exactly what the count-change bypasses (RealSeedCensusTests
/// established the seeds to N=11; this test does not re-run it, it exhibits the window it beats).</para>
///
/// <para>At each seed the coalescing pair is anchored to the census λ and its gap walked out along an offset
/// ladder; ONE artifact-free <see cref="EpCharacter"/> reading at the seed confirms the pair is a genuine
/// Jordan block. Grounds the F86a-vs-F89 reconciliation notes in PROOF_F86A_EP_MECHANISM / LocalGlobalEpLink.
/// Run: <c>dotnet test "compute/RCPsiSquared.Diagnostics.Tests" --filter "Category=F86A_MASKING"
/// --logger "console;verbosity=detailed"</c>.</para></summary>
public sealed class F86aSeedMaskingTests
{
    private readonly ITestOutputHelper _out;
    public F86aSeedMaskingTests(ITestOutputHelper o) => _out = o;

    // F86a's real scan: Q ∈ [0.5, 4.0], ~121 sweep points ⟹ ΔQ ≈ 0.0292, nearest grid point up to ≈0.0146 off.
    private const double F86aGridHalfStep = 0.0146;

    /// <summary>The gap between the two eigenvalues of the raw (1,2) block nearest the census λ (the coalescing
    /// √-EP pair, anchored by value so the ladder tracks the same pair), and their midpoint.</summary>
    private static (double Gap, Complex Center, Matrix<Complex> M) PairAt(int n, double q, double lamA)
    {
        var m = Matrix<Complex>.Build.DenseOfArray(SectorBlock.Build(n, 1, 2, new Complex(q, 0)));
        var near = m.Evd().EigenValues.ToArray()
                    .OrderBy(z => (z - new Complex(lamA, 0)).Magnitude).Take(2).ToArray();
        return ((near[0] - near[1]).Magnitude, 0.5 * (near[0] + near[1]), m);
    }

    private bool ReportSeed(int n, RealSeed s)
    {
        _out.WriteLine($"  seed q*≈{s.QStar:F6} (R{s.RParity:+0;-0}), census λ≈{s.LambdaA:F4}:");

        // ONE artifact-free character reading at the seed: is the coalescing pair a Jordan block?
        var (gap0, center, m) = PairAt(n, s.QStar, s.LambdaA);
        var spec = m.Evd().EigenValues.ToArray();
        double third = spec.Select(z => (z - center).Magnitude).OrderBy(d => d).Skip(2).First();
        double radius = Math.Min(0.45 * third, Math.Max(5.0 * gap0, 0.02));
        var reading = EpCharacter.Characterize(m, center, radius, quadPoints: 200);
        _out.WriteLine($"      at the seed: EpCharacter = {reading.Kind} (alg={reading.Algebraic}, geo={reading.Geometric}, " +
                       $"relDep={reading.Departure / Math.Max(1.0, reading.CompressionNorm):E1}), pair gap={gap0:E1}");

        // The window: the same pair's gap as we step away from q*.
        _out.WriteLine("      coalescing-pair gap vs distance from the seed:");
        double gapAtGrid = 0;
        foreach (double off in new[] { 0.0, 1e-3, 5e-3, F86aGridHalfStep })
        {
            double g = PairAt(n, s.QStar + off, s.LambdaA).Gap;
            _out.WriteLine($"        |q-q*|={off:F4}   gap={g:E2}" +
                           (Math.Abs(off - F86aGridHalfStep) < 1e-9 ? "   <- F86a grid half-step" : ""));
            if (Math.Abs(off - F86aGridHalfStep) < 1e-9) gapAtGrid = g;
        }
        bool defectiveHere = reading.Kind == EpCharacter.EpKind.Defective;
        bool wideAtGrid = gapAtGrid > 30 * gap0;   // the pair is plainly split one grid-step away
        _out.WriteLine($"      => Jordan block at the seed: {defectiveHere}; pair {gapAtGrid / Math.Max(gap0, 1e-30):F0}x " +
                       "wider at F86a's grid offset (unresolvable as a coalescence there).");
        return defectiveHere && wideAtGrid;
    }

    private int ReportN(int n)
    {
        _out.WriteLine($"=== N={n}, raw (1,2) block, dim {n * (n * (n - 1) / 2)} ===");
        // One clean census seed per N (the N=9 one is the SeedHolonomyWitness locus).
        double target = n == 5 ? 1.077615 : 0.849011;
        var s = RealDefectiveSeeds.ForN(n).OrderBy(x => Math.Abs(x.QStar - target)).First();
        int shown = ReportSeed(n, s) ? 1 : 0;
        _out.WriteLine("");
        return shown;
    }

    [Fact]
    [Trait("Category", "F86A_MASKING")]
    public void N5_And_N9_TheWindowThatHidTheSeedFromF86a()
    {
        int five = ReportN(5);
        int nine = ReportN(9);

        _out.WriteLine("Detection window ~1e-3 vs F86a grid step ~0.029 (a ~20-30x mismatch): a sit-and-characterize");
        _out.WriteLine("scan on that grid cannot land on the cusp, so F86a read 'no real defective EP'. F89's");
        _out.WriteLine("count-change detector is grid-robust (an integer real-root jump across q*), flags the");
        _out.WriteLine("transition at ΔQ=0.01 and bisects onto q*. Same block, same axis, a better instrument.");

        Assert.True(five > 0, "N=5: the seed must be a Jordan block AND plainly split one grid-step away.");
        Assert.True(nine > 0, "N=9: the seed must be a Jordan block AND plainly split one grid-step away.");
    }
}
